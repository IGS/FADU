#!/usr/bin/env julia

"""
FADU.jl - Feature Aggregate Depth Utility.

Description - Generate fractional counts of alignments that map to non-overlapping portions of genes
Requires (version number listed is earliest supported version):
    Julia - v0.7
    GenomicFeatures.jl - v1.0.0
    BioAlignments.jl - v1.0.0

By: Shaun Adkins (sadkins@som.umaryland.edu)
    Matthew Chung (mattchung@umaryland.edu)

"""

# The macro on modules and functions makes the code available to all worker processes
using ArgParse
using BioAlignments: BAM, SAM
using GenomicFeatures
using Printf

include("alignment_overlaps.jl")
include("bam_record.jl")
include("feature_counts.jl")
include("gff_feature.jl")

const VERSION_NUMBER = "1.6"    # Version number of the FADU program
const MAX_FRAGMENT_SIZE = 1000 # Maximum size of fragment.  If exceeded, fragment will be considered two reads
const EM_ITER_DEFAULT = 1 # Number of iterations to do EM-algorithm

### Functions to clean up, then move to their proper "include" script


########
# Main #
########

function parse_commandline()
    s = ArgParseSettings(description = "Generate counts of reads that map to non-overlapping portions of genes",
    prog = "fadu.jl",
	version = VERSION_NUMBER,
	add_version = true,
    add_help = true)

    # The macro to add a table of arguments and options to the ArgParseSettings object
    @add_arg_table s begin
        "--bam_file", "-b"
            help = "Path to BAM file (SAM is not supported)."
            metavar = "/path/to/file.bam"
            required = true
        "--gff3_file", "-g"
            help = "Path to GFF3-formatted annotation file (GTF is not supported)."
            metavar = "/path/to/annotation.gff3"
            required = true
        "--output_dir", "-o"
            help = "Directory to write the output."
            metavar = "/path/to/output/dir/"
            required = true
        "--stranded", "-s"
            help = "Indicate if BAM reads are from a strand-specific assay. Choose between 'yes', 'no', or 'reverse'."
            default = "no"
        "--feature_type", "-f"
            help = "Which GFF3/GTF feature type (column 3) to obtain. readcount statistics for. Case-sensitive."
            default = "gene"
        "--attribute_type", "-a"
            help = "Which GFF3/GTF feature type (column 9) to obtain. readcount statistics for. Case-sensitive."
            default = "ID"
        "--keep_only_proper_pairs", "-p"
            help = "If enabled, keep only properly paired reads when performing calculations."
            action = :store_true
            dest_name = "pp_only"
        "--max_fragment_size", "-m"
            help = "If the fragment size of properly-paired reads exceeds this value, process pair as single reads instead of as a fragment. Setting this value to 0 will make every fragment pair be processed as two individual reads. If --keep_only_proper_pairs is enabled, then any fragment exceeding this value will be discarded."
            default = MAX_FRAGMENT_SIZE
            arg_type = Int
            range_tester = (x->typemin(UInt)<=x<=typemax(UInt))
        "--remove_multimapped", "-M"
            help = "If enabled, remove any reads or fragments that are mapped to multiple regions of the genome, indiated by the 'NH' attribute being greater than 1.  Otherwise, use EM algorithm to quantify reads after all other reads are counted."
            action = :store_true
            dest_name = "rm_multimap"
        "--em_iterations", "-e"
            help = "Number of iterations to perform EM algorithm on multimapped read depth. Only applies if --remove_multimapped flag is not passed (is disabled)"
            default = EM_ITER_DEFAULT
            arg_type = Int
            range_tester = (x->x>0)
            dest_name = "em_iter"

    # Will not add log_file or debug options for now
    end
    # Converts the ArgParseSettings object into key/value pairs
    return parse_args(s)
end

function validate_args(args)
    """Validate the passed arguments."""
    isfile(args["gff3_file"]) || throw(FileNotFoundException("GFF3 file does not seem to exist. Please check supplied path."))
    isfile(args["bam_file"]) || throw(FileNotFoundException("BAM file does not seem to exist. Please check supplied path."))
    if !isdir(args["output_dir"])
        @debug("Creating output directory")
        mkdir(args["output_dir"])
    end
    args["stranded"] in ["yes", "no", "reverse"] || error("--stranded argument must be either 'yes', 'no', or 'reverse'.")
end

function main()
    args = parse_commandline()
    validate_args(args)
    @info("Parsed args:")
    for (arg,val) in args
        @info("  $arg  =>  $val")
    end

    @info("Processing annotation features...")
    features = open(collect, GFF3.Reader, args["gff3_file"])    # Array of GFF3.Record objects
    # Only keep features for the chosen feature type
    filter!(x -> GFF3.featuretype(x) == args["feature_type"], features)

    @info("Getting unique coordinates per contig and feature...")
    uniq_coords = create_uniq_coords_dict(features, args["stranded"])

    @info("Initializing dictionary of feature count information")
    feat_overlaps = create_feat_overlaps_dict(features, uniq_coords, args["attribute_type"], args["stranded"])

    bai_file = replace(args["bam_file"], ".bam" => ".bai")
    # Attempt to find index file for BAM file in same directory
    if !isfile(bai_file)
        bai_file = args["bam_file"] * ".bai"
        if !isfile(bai_file)
            bam = args["bam_file"]
            error("Attempted to find index file for BAM file $bam but could not find one")
        end
    end
    @debug("Found BAM index file at $bai_file")

    multimapped_dict = Dict{String, IntervalCollection}()

    # Open a BAM file and iterate over records overlapping mRNA transcripts.
    @info("Opening BAM alignment file...")
    reader = open(BAM.Reader, args["bam_file"], index = bai_file)
    @info("Now finding overlaps between alignment and annotation records...")
    for feature in features
        featurename = get_featurename_from_attrs(feature, args["attribute_type"])
        process_overlaps!(feat_overlaps[featurename], multimapped_dict, reader, feature, args)
    end
    close(reader)

    # If multimapped reads are kept, use EM algorithm to count and re-add back into the feature counts
    if !args["rm_multimap"]
        featurenames = collect(keys(feat_overlaps))
        mm_feat_overlaps = Dict{String, FeatureOverlap}(featurename => initialize_overlap_info(feat_overlaps[featurename].coords_set) for featurename in featurenames)

        @debug("Multimapped alignment templates: ", length(multimapped_dict))
        @info("Counting and adjusting multimapped alignment feature counts via Expectation-Maximization algorithm...")
        while args["em_iter"] > 0
            @debug("\tEM iterations left: ", args["em_iter"])
            adjusted_overlaps = merge_mm_counts(feat_overlaps, mm_feat_overlaps, false)
            args["em_iter"] -= 1
            mm_feat_overlaps = compute_mm_counts_by_em(adjusted_overlaps, multimapped_dict, features, args)
        end
        # Last iteration the alignment counts are added too
        merge_mm_counts!(feat_overlaps, mm_feat_overlaps, true)
    end

    # Calculate sum of all counts for each feature
    totalcounts = calc_totalcounts(feat_overlaps)

    @info("Writing counts output to file...")
    out_file = joinpath(args["output_dir"], splitext(basename(args["bam_file"]))[1]) * ".counts.txt"
    out_f = open(out_file, "w")
    write(out_f, "featureID\tuniq_len\tnum_alignments\tcounts\ttpm\n")
    # Write output
    for featurename in sort(collect(keys(feat_overlaps)))
        uniq_len::UInt = length(feat_overlaps[featurename].coords_set)
        num_alignments = feat_overlaps[featurename].num_alignments
        feat_counts = feat_overlaps[featurename].feat_counts
        if isnan(feat_counts)
            feat_counts = zero(Float32)
        end
        tpm = calc_tpm(uniq_len, totalcounts, feat_counts)
        if isnan(tpm)
            tpm = zero(Float32)
        end
        s = @sprintf("%s\t%i\t%.1f\t%.2f\t%.2f\n", featurename, uniq_len, num_alignments, feat_counts, tpm)
        write(out_f, s)
    end
    close(out_f)

    @info("FADU is complete!  Exiting!")
end

main()
exit()
