#!/usr/bin/env julia

"""
FADU.jl - Feature Aggregate Depth Utility.

Description - Generate fractional counts of alignments that map to non-overlapping portions of genes
Requires (version number listed is earliest supported version):
    Julia - v1.7
    * Packages *
    BGZFStreams.jl - v0.3.2
    GenomicFeatures.jl - v2.1.0
    GFF3 - v0.2.3
    Indexes.jl - v0.1.3
    StructArrays.jl - v0.6.17
    XAM.jl - v0.3.1
    BED.jl - v0.3.0

By: Shaun Adkins (sadkins@som.umaryland.edu)
    Matthew Chung (mattchung@umaryland.edu)

"""

# The macro on modules and functions makes the code available to all worker processes
using ArgParse
using Logging
using XAM: BAM, SAM
using GenomicFeatures
using GFF3
using Indexes
using Printf
using StructArrays
using BED

include("alignment_overlaps.jl")
include("bam_record.jl")
include("feature_counts.jl")
include("gff_feature.jl")

const VERSION_NUMBER = "1.9"    # Version number of the FADU program
const MAX_FRAGMENT_SIZE = 1000 # Maximum size of fragment.  If exceeded, fragment will be considered two reads
const EM_ITER_DEFAULT = 1 # Number of iterations to do EM-algorithm

### Functions to clean up, then move to their proper "include" script

function remove_excluded_regions!(uniq_coords::Dict{String, Dict}, excluded_regions_file, stranded::Bool=false)
    """Remove any coordinates that overlap the regions in the BED file. Edits "uniq_coords" in-place."""

    @info("Removing alignments that overlap excluded regions...")
    # if file does not exist, throw error
    isfile(excluded_regions_file) || throw(SystemError("Excluded regions file does not seem to exist. Please check supplied path."))

    exclude_regions = open(collect, BED.Reader, excluded_regions_file)
    exclude_coords = Dict{String, Dict}()
    # Read in all excluded regions, by strand if the "stranded" argument is passed
    for region in exclude_regions
        # Store chrom since that should match seqid
        seqid = BED.chrom(region)
        if !haskey(uniq_coords, seqid)
            continue
        end

        # "get!" will create the key if it doesn't exist in addition to returning the value
        seqid_exclude_coords = get!(exclude_coords, seqid, Dict{Char, BitSet}())

        start = BED.chromstart(region)
        stop = BED.chromend(region)

        # Store coordinates to exclude
        # If unstranded, then everything goes in "+"
        if stranded
            strand = convert(Char, BED.strand(region))
            union!(get!(seqid_exclude_coords, strand, BitSet()), BitSet(start:stop))
        else
            union!(get!(seqid_exclude_coords, "+", BitSet()), BitSet(start:stop))
        end
    end


    # Remove any alignments that overlap these regions
    for (seqid, strand_coords) in exclude_coords
        if haskey(uniq_coords, seqid)
            for (strand, coords) in strand_coords
                if haskey(uniq_coords[seqid], strand)
                    setdiff!(uniq_coords[seqid][strand], coords)
                end
            end
        end
    end

end


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
    @add_arg_table! s begin
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
        "--no_output_header"
            help = "If enabled, do not write the header line to the output file."
            action = :store_true
            dest_name = "no_header"
        "--exclude_regions_file", "-x"
            help = "Path to BED file containing regions to exclude from analysis.  Any alignments that overlap these regions will be ignored."
            metavar = "/path/to/regions.bed"
        "--stranded", "-s"
            help = "Indicate if BAM reads are from a strand-specific assay. Choose between 'yes', 'no', or 'reverse'."
            default = "no"
            range_tester = (x->x in ["yes", "no", "reverse"])
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
        "--log_level"
            help = "Set the log level.  Options are: DEBUG, INFO, WARNING, ERROR, CRITICAL."
            default = "INFO"
            metavar = "LEVEL"
            range_tester = (x->x in ["DEBUG", "INFO", "WARNING", "ERROR"])
        "--log_file"
            help = "Path to log file.  If not specified, log messages will be printed to STDERR."
            metavar = "/path/to/logfile.log"

    # Will not add log_file or debug options for now
    end
    # Converts the ArgParseSettings object into key/value pairs
    return parse_args(s)
end

function setup_logger(args)
    """Set up logging."""
    stream = args["log_file"] === nothing ? stderr : open(args["log_file"], "w")
    logleveldict = Dict("DEBUG" => -1000, "INFO" => 0, "WARNING" => 1000, "ERROR" => 2000)
    loglevel = get(logleveldict, uppercase(args["log_level"]), 0)
    logger = ConsoleLogger(stream, loglevel)
    global_logger(logger)
    # TODO: prevent debug messages from printing the file and line number

end

function validate_args(args)
    """Validate the passed arguments."""
    @info("Validating arguments...")
    isfile(args["bam_file"]) || error("BAM file does not exist. Please check supplied path")
    isfile(args["gff3_file"]) || error("GFF3 file does not exist. Please check supplied path")
    if !isdir(args["output_dir"])
        @info("Creating output directory at ", args["output_dir"])
        mkdir(args["output_dir"])
    end
end

function main()
    args = parse_commandline()
    setup_logger(args)

    try
        validate_args(args)
    catch e
        @error(e)
        exit(1)
    end

    @info("Parsed args:")
    for (arg,val) in args
        @info("  $arg  =>  $val")
    end

    @info("Processing annotation features...")
    features = try
        open(collect, GFF3.Reader, args["gff3_file"])    # Array of GFF3.Record objects
    catch e
        println(stderr, "ERROR: GFF3 file was unable to be read, due to possible improper formatting.  Please consult https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md for how to properly format the GFF3 file.  Actual error is below.  Program will now exit.\n")
        showerror(stderr, e)
        exit(1)
    end

    # Only keep features for the chosen feature type
    filter!(x -> GFF3.featuretype(x) == args["feature_type"], features)

    @info("Getting unique coordinates per contig and feature...")
    uniq_coords = create_uniq_coords_dict(features, args["stranded"])

    # If exclude_regions is passed, then remove any alignments that overlap these regions
    if args["exclude_regions_file"] !== nothing
        try
            remove_excluded_regions!(uniq_coords, args["exclude_regions_file"], isstranded(args["stranded"]))
        catch e
            println(stderr, "ERROR: Excluded regions file was unable to be read, due to possible improper formatting.  Please consult https://genome.ucsc.edu/FAQ/FAQformat.html#format1 for how to properly format the BED file.  Actual error is below.  Program will now exit.\n")
            showerror(stderr, e)
            exit(1)
        end
    end

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

    # Open a BAM file and iterate over records overlapping GFF features.
    @info("Opening BAM alignment file...")
    reader = open(BAM.Reader, args["bam_file"], index = bai_file)

    @info("Now finding overlaps between alignment and annotation records...")
    multimapped_dict = process_all_feature_overlaps(feat_overlaps, features, reader, args)
    close(reader)

    # If multimapped reads are kept, use EM algorithm to count and re-add back into the feature counts
    if !args["rm_multimap"]
        num_multimaps = length(collect(keys(multimapped_dict)))
        @debug("Multimapped alignment templates: ", num_multimaps)

        featurenames = collect(keys(feat_overlaps))
        mm_feat_overlaps = Dict{String, FeatureOverlap}(featurename => initialize_overlap_info(coordinate_set(feat_overlaps[featurename])) for featurename in featurenames)

        @info("Counting and adjusting multimapped alignment feature counts via Expectation-Maximization algorithm...")
        while args["em_iter"] > 0
            @debug("\tEM iterations left: ", args["em_iter"])
            adjusted_overlaps = merge_mm_counts(feat_overlaps, mm_feat_overlaps, false)
            args["em_iter"] -= 1
            mm_feat_overlaps = compute_mm_counts_by_em(adjusted_overlaps, multimapped_dict)
        end
        # Last iteration the alignment counts are added too
        merge_mm_counts!(feat_overlaps, mm_feat_overlaps, true)
    end

    # Calculate sum of all counts for each feature
    totalcounts = calc_totalcounts(feat_overlaps)

    @info("Writing counts output to file...")
    out_file = joinpath(args["output_dir"], splitext(basename(args["bam_file"]))[1]) * ".counts.txt"
    out_f = open(out_file, "w")
    if !args["no_header"]
        write(out_f, "featureID\tuniq_len\tnum_alignments\tcounts\ttpm\n")
    end
    # Write output, sorted alphabetically
    for featurename in sort(collect(keys(feat_overlaps)))
        uniq_len::UInt = length(coordinate_set(feat_overlaps[featurename]))
        num_alignments = totalalignments(feat_overlaps[featurename])
        feat_counts = featurecounts(feat_overlaps[featurename])
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
