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

const VERSION_NUMBER = "1.3"    # Version number of the FADU program
const MAX_FRAGMENT_SIZE = 1000 # Maximum size of fragment.  If exceeded, fragment will be considered two reads
const MIN_MAP_QUAL = 10 # Minimum mapping quality
const CHUNK_SIZE = 10000000 # Number of valid BAM fragments to read in before determining overlaps
# NOTE: Making chunk_counter a UInt32, so this constant should not exceed 4,294,967,295

#is_duplicate(record::BAM.Record) = BAM.flag(record) & SAM.FLAG_DUP == 0x0400
is_mate_reverse(record::BAM.Record) = BAM.flag(record) & SAM.FLAG_MREVERSE == 0x0020
is_proper_pair(record::BAM.Record) = BAM.flag(record) & SAM.FLAG_PROPER_PAIR == 0x0002
is_read1(record::BAM.Record) = BAM.flag(record) & SAM.FLAG_READ1 == 0x0040
is_read2(record::BAM.Record) = BAM.flag(record) & SAM.FLAG_READ2 == 0x0080
is_reverse(record::BAM.Record) = BAM.flag(record) & SAM.FLAG_REVERSE == 0x0010

# Forward-stranded assay
## Positive Strand:
### R1 - forward (65), R2 - revcom (145)
## Negative Strand
### R1 - revcom (81), R2 - forward (129)
# Reverse-stranded assay (i.e. Illumina) have these strands flipped

### For non-properly paired reads
# Singletons will either have flag 0 or 16

# Forward strand flags
flag_65(record::BAM.Record) = is_read1(record) && !is_reverse(record)
flag_145(record::BAM.Record) = is_read2(record) && is_reverse(record)
# Reverse strand flags
flag_81(record::BAM.Record) = is_read1(record) && is_reverse(record)
flag_129(record::BAM.Record) = is_read2(record) && !is_reverse(record)


function add_nonoverlapping_feature_coords!(uniq_coords::Dict{String, Dict}, feature::Interval{GenomicFeatures.GFF3.Record}, stranded::Bool)
    """Adjust the set of coordinates stored in the Dict to be the symmetric difference with those in the feature."""
    seqid = seqname(feature)
    strand = get_strand_of_interval(feature, stranded)
    # First passthru of new contig or region, add first set of coords (since all uniq at this moment)
    if !(haskey(uniq_coords, seqid))
        uniq_coords[seqid] = Dict{Char,Set{UInt}}(strand => get_feature_coords_set(feature))
        return
    elseif !(haskey(uniq_coords[seqid], strand))
        uniq_coords[seqid][strand] = get_feature_coords_set(feature)
        return
    end
    curr_contig_coords = uniq_coords[seqid][strand]
    new_feat_coords = get_feature_coords_set(feature)
    # Get overlapping coordinates and remove from current contig coords set (symmetric difference)
    uniq_coords[seqid][strand] = symdiff(curr_contig_coords, new_feat_coords)
end

function assign_read_to_strand(record::BAM.Record, reverse_strand::Bool=false)
    """Use the bitwise flags to assign the paired read to the correct strand."""

    pos_flags = Set{Bool}([flag_145(record), flag_65(record)])
    neg_flags = Set{Bool}([flag_81(record), flag_129(record)])
    if reverse_strand
        pos_flags = Set{Bool}([flag_81(record), flag_129(record)])
        neg_flags = Set{Bool}([flag_145(record), flag_65(record)])
    end

    # Does record flag belong in the pos or neg set?
    any(pos_flags) && return '+'
    any(neg_flags) && return '-'

    # Read must be a singleton
    if is_reverse(record)
        # If reverse-stranded, strand is flipped
        reverse_strand && return '+'
        return '-'
    end
    reverse_strand && return '-'
    return '+'
end

function calc_tpm(len::UInt, depth_sum::Float32, feat_depth::Float32)
    """Calculate TPM score for current feature."""
    return @fastmath(feat_depth *1000 / len) * 1000000 / depth_sum
end

function compute_frag_feat_ratio(uniq_coords::Dict{String, Dict}, fragment::Interval{Char}, feature::GFF3.Record, strand::Char)
    """Calculate the ratio of fragament coordinates that intersect with non-overlapping feature coordinates."""
        # Pertinent fragment info
        frag_start = leftposition(fragment)
        frag_end = rightposition(fragment)
        frag_intersect = intersect(frag_start:frag_end, uniq_coords[GFF3.seqid(feature)][strand])
        feature_name = get_feature_name_from_attrs(feature, "ID")
        # Percentage of fragment that aligned with this annotation feature
        return @fastmath length(frag_intersect) / length(frag_start:frag_end)
end

function get_feature_coords_set(feature::Interval{GenomicFeatures.GFF3.Record})
    """Get the range of coordinates for this feature, returned as a Set."""
    return Set{UInt}(leftposition(feature) : rightposition(feature))
end

function get_feature_name_from_attrs(feature::GFF3.Record, attr_type::String)
    """Get attribute ID to use as the feature name."""
    gene_vector = GFF3.attributes(feature, attr_type)    # col 9
    validate_feature_attribute(gene_vector)
    return gene_vector[1]
end

function get_feature_nonoverlapping_length(feature::Interval{GenomicFeatures.GFF3.Record}, uniq_coords::Dict{String, Dict}, stranded::Bool)
    """Get number of nonoverlapping coordinates for given feature."""
    seqid = seqname(feature)
    strand = get_strand_of_interval(feature, stranded)
    feat_coords = get_feature_coords_set(feature)
    return length(intersect(feat_coords, uniq_coords[seqid][strand]))
end


function get_fragment_start_end(record::BAM.Record, record_type::Char)
    """Get the start and end coordinates of the fragment/read."""
    if record_type == 'R'
        return BAM.position(record):BAM.rightposition(record)
    end
    return BAM.nextposition(record):BAM.rightposition(record)
end

function get_fragment_interval(record::BAM.Record, record_type::Char, reverse_strand::Bool=false)
    """Return a fragment-based Interval for the current record.""" 
    return Interval(BAM.refname(record), get_fragment_start_end(record, record_type), assign_read_to_strand(record, reverse_strand), record_type)
end

function get_strand_of_interval(interval::Interval, stranded::Bool)
    """Get strand of interval with respect to strandedness arguments."""
    strand = '+'
    if stranded
        strand = convert(Char, interval.strand)
    end
    return strand
end

function increment_feature_overlap_information(feat_dict::Dict{String,Union{UInt,Float32}}, frag_feat_ratio::Float32, is_read::Bool)
    """Increment counter and depth information for feature if fragment overlapped with uniq coords."""
    if frag_feat_ratio > 0.0
        counter::Float32 = 1.0
        if is_read
            # If dealing with read, do not want to give as much weight since both reads will be added independently compared to a single fragment
            counter = 0.5
        end
        feat_dict["counter"] += counter
        feat_dict["feat_depth"] += frag_feat_ratio * counter
    end
end

function is_chunk_ready(counter::UInt32, chunk_size::UInt32)
    """Test to see if chunk is ready for further processing."""
    counter % chunk_size == 0 && return true
    return false
end

function is_multimapped(record::BAM.Record)
    """Test to see if alignment is multimapped across multiple regions of the genome."""
    try
        return haskey(record, "NH") && record["NH"] > 1
    catch
        # Ran into bug where some attributes are 'X' type which is not valid
        # Bioalignments BAM.Record.auxdata(record) throws LoadError for these
        return false
    end
end

function is_stranded(strand_type::String)
    """Check to see if stranded type is stranded or not."""
    strand_type == "no" && return false
    return true
end

function is_reverse_stranded(strand_type::String)
    """Check to see if stranded type is reverse-stranded."""
    strand_type == "reverse" && return true
    return false
end

function is_templength_negative(templength::Int64)
    """Check to see if the read template is going in the opposite version."""
    return templength < 0
end

function is_templength_smaller_than_max_fragment_size(templength::Int64, max_frag_size::UInt16)
    """Check to see if fragment template length is smaller than specified maximum fragment size."""
    return abs(templength) <= max_frag_size
end

function process_overlaps!(feat_overlaps::Dict{String, Dict}, uniq_coords::Dict{String, Dict}, fragment_intervals::IntervalCollection{Char}, features::IntervalCollection{GFF3.Record}, args::Dict)
    """Process current chunk of fragment intervals that overlap with feature intervals."""
    for (fragment, feature) in eachoverlap(fragment_intervals, features)
        feat_record = metadata(feature)
        GFF3.featuretype(feat_record) == args["feature_type"] || continue
        feature_name = get_feature_name_from_attrs(feat_record, args["attribute_type"])

        bam_strand = get_strand_of_interval(fragment, is_stranded(args["stranded"]))
        gff_strand = get_strand_of_interval(feature, is_stranded(args["stranded"]))
        if bam_strand == gff_strand
            frag_feat_ratio::Float32 = compute_frag_feat_ratio(uniq_coords, fragment, feat_record, bam_strand)
            increment_feature_overlap_information(feat_overlaps[feature_name], frag_feat_ratio, metadata(fragment)=='R')
        end
    end
end

function validate_feature_attribute(gene_vector::Vector{String})
    """Ensure only one attribute of this type exists for this feature record."""
    length(gene_vector) > 0 || error("ERROR - Attribute field 'ID' found to have no entries.\n")
    length(gene_vector) == 1 || error("ERROR - Attribute field 'ID' found to have multiple entries.\n")
end

function validate_fragment(record::BAM.Record)
    """Ensure alignment record can be used to calculate fragment depth. (Only need negative template of fragment)"""
    return is_proper_pair(record) && is_templength_negative(BAM.templength(record))
end

function validate_read(record::BAM.Record, max_frag_size::UInt16)
    """Ensure alignment record can be used to calculate read depth."""
    return !(is_proper_pair(record) && is_templength_smaller_than_max_fragment_size(BAM.templength(record), max_frag_size))
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
            help = "If the fragment size of properly-paired reads exceeds this value, process pair as single reads instead of as a fragment. Setting this value to 0 will make every fragment pair be processed as two individual reads. Maximum value is 65535 (to allow for use of UInt16 type, and fragments typically are not that large). If --keep_only_proper_pairs is enabled, then any fragment exceeding this value will be discarded."
            default = MAX_FRAGMENT_SIZE
            arg_type = Int
            range_tester = (x->typemin(UInt16)<=x<=typemax(UInt16))
        "--min_mapping_quality", "-q"
            help = "Set the minimum mapping quality score for a read to be considered. Max value accepted is 255."
            default = MIN_MAP_QUAL
            arg_type = Int
            range_tester = (x->typemin(UInt8)<=x<=typemax(UInt8))
        "--remove_multimapped", "-M"
            help = "If enabled, remove any reads or fragments that are mapped to multiple regions of the genome, indiated by the 'NH' attribute being greater than 1."
            action = :store_true
            dest_name = "rm_multimap"
        "--chunk_size", "-C"
            help = "Number of validated reads to store into memory before processing overlaps with features."
            default = CHUNK_SIZE
            arg_type = Int
            range_tester = (x->typemin(UInt32)<=x<=typemax(UInt32))

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
        @debug "Creating output directory" 
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
    gff3_reader = open(GFF3.Reader, args["gff3_file"])
    features = IntervalCollection(gff3_reader)
    close(gff3_reader)

    @info("Getting unique coordinates per contig and feature...")
    uniq_coords = Dict{String, Dict}()
    for feature in features
        GFF3.featuretype(metadata(feature)) == args["feature_type"] || continue
        add_nonoverlapping_feature_coords!(uniq_coords, feature, is_stranded(args["stranded"]))
    end
    @info("Initializing feature overlap dictionary...")
    feat_overlaps = Dict{String, Dict}()
    for feature in features
        GFF3.featuretype(metadata(feature)) == args["feature_type"] || continue
        feature_name = get_feature_name_from_attrs(metadata(feature), args["attribute_type"])
        feature_len::UInt = get_feature_nonoverlapping_length(feature, uniq_coords, is_stranded(args["stranded"]))
        counter::Float32 = 0.0; feat_depth::Float32 = 0.0
        feat_overlaps[feature_name] = Dict{String, Union{UInt,Float32}}("counter" => counter, "feat_depth" => counter, "uniq_len" => feature_len)
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
    record = BAM.Record()
    fragment_intervals = IntervalCollection{Char}()
    valid_record_counter::UInt32 = 0
    max_frag_size = convert(UInt16, args["max_fragment_size"])
    chunk_size = convert(UInt32, args["chunk_size"])
    record_type = '0'
    while !eof(bam_reader)
        read!(bam_reader, record)
        # Validation of current read
        BAM.ismapped(record) || continue
        BAM.isprimary(record) || continue
        BAM.mappingquality(record) >= args["min_mapping_quality"] || continue 
        args["rm_multimap"] && is_multimapped(record) && continue
        if validate_fragment(record) && is_templength_smaller_than_max_fragment_size(BAM.templength(record), max_frag_size)
            record_type = 'F'
        elseif !args["pp_only"] && validate_read(record, max_frag_size)
            record_type = 'R'
        else
            continue
        end
        if record_type == '0'
            error("Record passed through validation without being assigned 'R' or 'F'\n $record")
        end

        # Store reads as chunks and process chunks when chunk is max size
        valid_record_counter += 1
        if is_chunk_ready(valid_record_counter, chunk_size)
            process_overlaps!(feat_overlaps, uniq_coords, fragment_intervals, features, args)
            fragment_intervals = IntervalCollection{Char}()
        end
        push!(fragment_intervals, get_fragment_interval(record, record_type, is_reverse_stranded(args["stranded"])))
    end
    # Final chunk
    process_overlaps!(feat_overlaps, uniq_coords, fragment_intervals, features, args)
    close(bam_reader)

    # Calculate sum of all counts for each feature
    depth_sum = sum(feat_overlaps[feat]["feat_depth"] for feat in keys(feat_overlaps))

    @info("Writing counts output to file...")
    out_file = joinpath(args["output_dir"], splitext(basename(args["bam_file"]))[1]) * ".counts.txt"
    out_f = open(out_file, "w")
    write(out_f, "featureID\tuniq_len\tnum_alignments\tcounts\ttpm\n")
    # Write output
    for feat_id in sort(collect(keys(feat_overlaps)))
        uniq_len = feat_overlaps[feat_id]["uniq_len"]
        counter = feat_overlaps[feat_id]["counter"]
        feat_depth = feat_overlaps[feat_id]["feat_depth"]
        tpm = calc_tpm(uniq_len, depth_sum, feat_depth)
        s = @sprintf("%s\t%i\t%.1f\t%.2f\t%.2f\n", feat_id, uniq_len, counter, feat_depth, tpm)
        write(out_f, s)
    end
    close(out_f)

    @info("FADU is complete!  Exiting!")
end

main()
exit()
