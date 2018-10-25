#!/usr/bin/env julia

"""
FADU.jl - Feature Aggregate Depth Utility.

Description - Generate fractional counts of fragments that map to non-overlapping portions of genes
Requires (version number listed is earliest supported version): 
    Julia - v0.7
    GenomicFeatures.jl - v1.0.0
    BioAlignments.jl - v1.0.0

By: Shaun Adkins (sadkins@som.umaryland.edu)
    Matthew Chung (mattchung@umaryland.edu)

"""

using ArgParse
using BioAlignments
using GenomicFeatures
using Printf

const VERSION_NUMBER = "1.0"    # Version number of the FADU program
const CHUNK_SIZE = 10000000 # Number of valid BAM fragments to read in before determining overlaps

#is_duplicate(record::BAM.Record) = BAM.flag(record) & SAM.FLAG_DUP == 0x0400
is_mate_reverse(record::BAM.Record) = BAM.flag(record) & SAM.FLAG_MREVERSE == 0x0020
is_proper_pair(record::BAM.Record) = BAM.flag(record) & SAM.FLAG_PROPER_PAIR == 0x0002
is_read1(record::BAM.Record) = BAM.flag(record) & SAM.FLAG_READ1 == 0x0040
#is_read2(record::BAM.Record) = BAM.flag(record) & SAM.FLAG_READ2 == 0x0080
is_reverse(record::BAM.Record) = BAM.flag(record) & SAM.FLAG_REVERSE == 0x0010

function add_nonoverlapping_feature_coords!(uniq_coords::Dict{String, Dict}, feature::Interval{GenomicFeatures.GFF3.Record}, stranded::Bool)
    """Adjust the set of coordinates stored in the Dict to be the symmetric difference with those in the feature."""
    seqid = seqname(feature)
    strand = '+'
    if stranded
        strand = convert(Char, feature.strand)
    end
    # First passthru of new contig or region, add first set of coords (since all uniq at this moment)
    if !(haskey(uniq_coords, seqid))
        uniq_coords[seqid] = Dict{Char,Set{Int}}(strand => Set{Int}(leftposition(feature) : rightposition(feature)))
        return
    elseif !(haskey(uniq_coords[seqid], strand))
        uniq_coords[seqid][strand] = Set{Int}(leftposition(feature) : rightposition(feature))
        return
    end
    curr_contig_coords = uniq_coords[seqid][strand]
    new_feat_coords = Set{Int}(leftposition(feature) : rightposition(feature))
    # Get overlapping coordinates and remove from current contig coords set (symmetric difference)
    uniq_coords[seqid][strand] = symdiff(curr_contig_coords, new_feat_coords)
end

function assign_read_to_strand(record::BAM.Record, reverse_strand=false)
    """Use the bitwise flags to assign the paired read to the correct strand."""
    # Forward-stranded assay
    ## Positive Strand:
    ### R1 - forward (97), R2 - revcom (145)
    ## Negative Strand
    ### R1 - revcom (81), R2 - forward (161)
    # Reverse-stranded assay (i.e. Illumina)
    ## Positive Strand
    ### R1 - revcom (81), R2 - forward (161)
    ## Negative Strand
    ### R1 - forward (97), R2 - revcom (145)

    # Forward strand flags (ignoring read2 flags)
    flag_97 = is_read1(record) && is_mate_reverse(record)
    flag_65 = is_read1(record) && !(is_reverse(record) || is_mate_reverse(record))
    # Reverse strand flags
    flag_81 = is_read1(record) && is_reverse(record)
    flag_113 = is_read1(record) && is_reverse(record) && is_mate_reverse(record)

    pos_flags = Set([flag_97, flag_65])
    neg_flags = Set([flag_81, flag_113])
    if reverse_strand
        pos_flags = Set([flag_81, flag_113])
        neg_flags = Set([flag_97, flag_65])
    end

    # Does record flag belong in the pos or neg set?
    any(pos_flags) && return '+'
    any(neg_flags) && return '-'

    # Most likely will never reach this
    query_name = BAM.tempname(record)
    @warn("Read $query_name did not get assigned to either strand apparently.")
    return '?'
end

function compute_frag_feat_ratio(uniq_coords::Dict{String, Dict}, fragment::Interval{String}, feature::GFF3.Record, stranded::Bool)
    """Calculate the ratio of fragament coordinates that intersect with non-overlapping feature coordinates."""
        # Pertinent fragment info
        strand = '+'
        if stranded
            strand = convert(Char, GFF3.strand(feature))
        end
        frag_start = leftposition(fragment)
        frag_end = rightposition(fragment)
        frag_intersect = intersect(frag_start:frag_end, uniq_coords[GFF3.seqid(feature)][strand])
        # Percentage of fragment that aligned with this annotation feature
        return length(frag_intersect) / length(frag_start:frag_end)
end

function get_feature_name_from_attrs(feature::GFF3.Record, attr_type::String)
    """Get attribute ID to use as the feature name."""
    gene_vector = GFF3.attributes(feature, attr_type)    # col 9
    validate_feature_attribute(gene_vector)
    return gene_vector[1]
end

function get_fragment_start_end(record::BAM.Record)
    """Get the start and end coordinates of the fragment."""
    frag_start = BAM.position(record)
    frag_end = frag_start + BAM.templength(record) - 1
    # If read is reverse... 
    if BAM.templength(record) < 0
        frag_end = BAM.rightposition(record)
        frag_start = frag_end + BAM.templength(record) + 1
    end
    return frag_start:frag_end
end

function get_fragment_interval(record::BAM.Record, reverse_strand::Bool=false)
    """Return a fragment-based Interval for the current record.""" 
    return Interval(BAM.refname(record), get_fragment_start_end(record), assign_read_to_strand(record, reverse_strand), BAM.tempname(record))
end

function is_chunk_ready(counter::Int)
    counter % CHUNK_SIZE == 0 && return true
    return false
end

function is_stranded(strand_type::String)
    """Check to see if stranded type is reverse-stranded"""
    strand_type == "no" && return false
    return true
end

function is_reverse_stranded(strand_type::String)
    """Check to see if stranded type is reverse-stranded"""
    strand_type == "reverse" && return true
    return false
end

function process_overlaps!(feat_overlaps::Dict{String, Dict}, uniq_coords::Dict{String, Dict}, fragment_intervals::IntervalCollection{String}, features::IntervalCollection{GFF3.Record}, args::Dict)
    """Process current chunk of fragment intervals that overlap with feature intervals."""
    for (fragment, feature) in eachoverlap(fragment_intervals, features)
        # Pertinent feature info
        feat_record = metadata(feature)
        GFF3.featuretype(feat_record) == args["feature_type"] || continue
        feature_name = get_feature_name_from_attrs(feat_record, args["attribute_type"])
        # Initialize feature_name into overlaps Dict if key does not exist
        get!(feat_overlaps, feature_name, Dict{String, Real}("counter" => 0, "gene_depth" => 0))

        frag_feat_ratio = compute_frag_feat_ratio(uniq_coords, fragment, feat_record, is_stranded(args["stranded"]))
        if frag_feat_ratio > 0
            feat_overlaps[feature_name]["counter"] += 1
            feat_overlaps[feature_name]["gene_depth"] += frag_feat_ratio
        end
    end
end

function validate_feature_attribute(gene_vector::Vector{String})
    """Ensure only one attribute of this type exists for this feature record."""
    length(gene_vector) > 0 || error("ERROR - Attribute field 'ID' found to have no entries.\n")
    length(gene_vector) == 1 || error("ERROR - Attribute field 'ID' found to have multiple entries.\n")
end

function validate_record(record::BAM.Record)
    """Ensure alignment record is a primary alignment, is not multimapped, and is properly paired."""
    return BAM.isprimary(record) && is_proper_pair(record)
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
    @add_arg_table s begin
        "--bam_file", "-b"
            help = "Path to BAM file (SAM is not supported)."
            metavar = "/path/to/file.bam"
            required = true
            #group = "bam group"
        "--gff3_file", "-g"
            help = "Path to GFF3-formatted annotation file (GTF is not supported)."
            metavar = "/path/to/annotation.gff3"
            required = true
            #group = "annot group"
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
        #"--count_by", "-c"
        #    help = "How to count the reads when performing depth calculations."
        #    default = "read"
        # "--keep_only_properly_paired"
        #     help = "Enable flag to remove any reads that are not properly paired from the depth count statistics"
        #     action = :store_true
        # "--rm_multimapped_reads"
        #     help = "Enable flag to remove any multimapped reads ('NH' tag > 1)"
        #     action = :store_true

    # Will not add log_file or debug options for now
    end

    # Converts the ArgParseSettings object into key/value pairs
    return parse_args(s)
end

function validate_args(args)
    """Validate the passed arguments."""
    isfile(args["gff3_file"]) || throw(FileNotFoundException("GFF3 file does not seem to exist. Please check supplied path."))
    isfile(args["gff3_file"]) || throw(FileNotFoundException("BAM file does not seem to exist. Please check supplied path."))
    if !isdir(args["output_dir"])
        @debug "Creating output directory" 
        mkdir(args["output_dir"])
    end
    args["stranded"] in ["yes", "no", "reverse"] || error("--stranded argument must be either 'yes', 'no', or 'reverse'.")
    #args["count_by"] in ["read", "fragment"] || error("--count_by argument must be either 'read' or 'fragment'.")
end

function main()
    args = parse_commandline()
    validate_args(args)
    @info("Parsed args:")
    for (arg,val) in args
        @info("  $arg  =>  $val")
    end

    @info("Processing annotation features...")
    gff3_reader = open(GFF3.Reader, args["gff3_file"])
    @time features = IntervalCollection(gff3_reader)

    @info("Getting unique coordinates per contig...")
    uniq_coords = Dict{String, Dict}()
    @time for feature in features
        GFF3.featuretype(metadata(feature)) == args["feature_type"] || continue
        add_nonoverlapping_feature_coords!(uniq_coords, feature, is_stranded(args["stranded"]))
    end

    @info("Processing BAM alignment records...")
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
    bam_reader = open(BAM.Reader, args["bam_file"], index=bai_file)

    @info("Now finding overlaps between alignment and annotation records...")
    feat_overlaps = Dict{String, Dict}()

    record = BAM.Record()
    fragment_intervals = IntervalCollection{String}()
    valid_record_counter = 0
    @time while !eof(bam_reader)
        read!(bam_reader, record)
        validate_record(record) && is_read1(record) || continue
        valid_record_counter += 1
        if is_chunk_ready(valid_record_counter)
            process_overlaps!(feat_overlaps, uniq_coords, fragment_intervals, features, args)
            fragment_intervals = IntervalCollection{String}()
        end
        push!(fragment_intervals, get_fragment_interval(record, is_reverse_stranded(args["stranded"])))
    end
    # Final chunk
    process_overlaps!(feat_overlaps, uniq_coords, fragment_intervals, features, args)
    close(bam_reader)

    @info("Determining which features did not have fragment alignments overlap...")
    lonely_features = []
    for record in gff3_reader
        GFF3.featuretype(record) == args["feature_type"] || continue
        feature_name = get_feature_name_from_attrs(record, args["attribute_type"])
        get!(feat_overlaps, feature_name, Dict{String, Real}("counter" => 0, "gene_depth" => 0))
    end

    @info("Writing counts output to file...")
    out_file = joinpath(args["output_dir"], splitext(basename(args["bam_file"]))[1]) * ".counts.txt"
    out_f = open(out_file, "w")
    write(out_f, "gene\tnum_alignments\tcounts\n")
    # Write gene counts and depth to file
    for gene_id in sort(collect(keys(feat_overlaps)))
        counter = feat_overlaps[gene_id]["counter"]
        gene_depth = feat_overlaps[gene_id]["gene_depth"]
        s = @sprintf("%s\t%i\t%.2f\n", gene_id, counter, gene_depth)
        write(out_f, s)
    end
    close(out_f)

    @info("FADU is complete!  Exiting!")
end

main()
exit()
