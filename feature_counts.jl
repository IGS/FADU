#!/usr/bin/env julia

"""
feature_counts.jl - Houses functions that deal with GFF feature-BAM alignment overlaps and counts

By: Shaun Adkins (sadkins@som.umaryland.edu)
"""

mutable struct FeatureOverlap
    num_alignments::Float32
    feat_counts::Float32
    coords_set::Set{UInt}
end

function adjust_mm_counts_by_em(mm_overlaps::Dict{String, FeatureOverlap}, feat_overlaps::Dict{String, FeatureOverlap}, alignment_dict::Dict{String,IntervalCollection}, features::Array{GFF3.Record,1}, args::Dict)
    """Adjust the feature counts of the multimapped overlaps via the Expectation-Maximization algorithm."""
    new_mm_overlaps = Dict{String, FeatureOverlap}()
    feat_and_mm_overlaps = copy(feat_overlaps)
    for feature_name in keys(feat_overlaps)
        new_mm_overlaps[feature_name] = initialize_overlap_info(feat_overlaps[feature_name].coords_set)
        new_mm_overlaps[feature_name].num_alignments = mm_overlaps[feature_name].num_alignments
        # Combine singly-mapped and multimapped feature counts
        feat_and_mm_overlaps[feature_name].feat_counts += mm_overlaps[feature_name].feat_counts
    end

    for record_tempname in keys(alignment_dict)
        # Get all features that align with the given multimapped template name on the same strand
        aln_feat_overlaps = Array{GenomicFeatures.GFF3.Record,1}()
        for aln_interval in alignment_dict[record_tempname]
            append!(aln_feat_overlaps, filter_alignment_feature_overlaps(features, aln_interval, isstranded(args["stranded"])))
        end

        template_totalcounts = calc_template_totalcounts(aln_feat_overlaps, feat_and_mm_overlaps, args["attribute_type"])
        
        # Adjust contribution proportion of multimapped reads, based on estimated relative abundance
        for aln_interval in alignment_dict[record_tempname]
            for feature in aln_feat_overlaps
                feature_name = get_feature_name_from_attrs(feature, args["attribute_type"])
                align_feat_ratio::Float32 = compute_align_feat_ratio(mm_overlaps[feature_name].coords_set, aln_interval)
                align_feat_ratio > 0.0 || continue
                aln_type = gettype_alignment(metadata(aln_interval))
                featurecount = align_feat_ratio * aln_type.count_multiplier
                feature_totalcounts = feat_and_mm_overlaps[feature_name].feat_counts
                rel_abundance = @fastmath feature_totalcounts / template_totalcounts
                new_mm_overlaps[feature_name].feat_counts += rel_abundance * featurecount
            end
        end
    end
    return new_mm_overlaps
end

function calc_template_totalcounts(aln_feat_overlaps::Array{GFF3.Record,1}, feat_and_mm_overlaps::Dict{String, FeatureOverlap}, attribute_type::String)
    """For all features overlapping the alignment template, return the sum of all those feature counts."""
    totalcounts = 0
    for feature in aln_feat_overlaps
        feature_name = get_feature_name_from_attrs(feature, attribute_type)
        totalcounts += feat_and_mm_overlaps[feature_name].feat_counts
    end
    return totalcounts
end

function calc_tpm(len::UInt, totalcounts::Float32, feat_counts::Float32)
    """Calculate TPM score for current feature."""
    return @fastmath(feat_counts *1000 / len) * 1000000 / totalcounts
end

function calc_totalcounts(feat_overlaps::Dict{String, FeatureOverlap})
    """Calculate the sum of all the feature counts."""
    return sum(feat_overlaps[feat].feat_counts for feat in keys(feat_overlaps))
end

function compute_align_feat_ratio(uniq_feat_coords::Set{UInt}, alignment::Interval)
    """Calculate the ratio of fragament coordinates that intersect with non-overlapping feature coordinates."""
        # Pertinent alignment info
        alignment_coords =  get_alignment_coords_set(alignment)
        alignment_intersect = intersect(uniq_feat_coords, alignment_coords)
        # Percentage of alignment that aligned with this annotation feature
        return @fastmath length(alignment_intersect) / length(alignment_coords)
end

function create_feat_overlaps_dict(features::Array{GenomicFeatures.GFF3.Record,1}, uniq_coords::Dict{String,Dict}, attr_type::String, stranded_type::String)
    feat_overlaps = Dict{String, FeatureOverlap}()
    for feature in features
        feature_name = get_feature_name_from_attrs(feature, attr_type)
        uniq_feat_coords = get_feature_nonoverlapping_coords(feature, uniq_coords, isstranded(stranded_type))
        feat_overlaps[feature_name] = initialize_overlap_info(uniq_feat_coords)
    end
    return feat_overlaps
end

function filter_alignment_feature_overlaps(features::Array{GenomicFeatures.GFF3.Record,1}, alignment::Interval, stranded::Bool)
    """Filter features to those just that align with the current alignment on the same strand."""
    aln_strand = getstrand(alignment, stranded)
    aln_feat_overlaps = filter(x -> isoverlapping(convert(Interval, x), alignment), features)
    filter!(x -> aln_strand == getstrand(x, stranded), aln_feat_overlaps)
    return aln_feat_overlaps
end

function increment_feature_overlap_information!(feat_overlap::FeatureOverlap, align_feat_ratio::Float32, aln_type::T) where {T<:AbstractAlignment}
    """Increment number of alignments and feature count information for feature if alignment overlapped with uniq coords."""
    feat_overlap.num_alignments += aln_type.count_multiplier
    feat_overlap.feat_counts += align_feat_ratio * aln_type.count_multiplier
end

function initialize_overlap_info(uniq_feat_coords::Set{UInt})
    """Initialize overlap diction information."""
    num_alignments::Float32 = 0.0; feat_counts::Float32 = 0.0
    return FeatureOverlap(num_alignments, feat_counts, uniq_feat_coords)
end

function merge_mm_counts!(feat_overlaps::Dict{String, FeatureOverlap}, mm_feat_overlaps::Dict{String, FeatureOverlap})
    """Merge EM-computed counts for multimapped reads back into the general feature overlap counts."""
    for feature in keys(feat_overlaps)
        feat_overlaps[feature].num_alignments += mm_feat_overlaps[feature].num_alignments
        feat_overlaps[feature].feat_counts += mm_feat_overlaps[feature].feat_counts
    end
end

function process_aln_interval_for_overlaps!(feat_overlaps::Dict{String, FeatureOverlap}, features::Array{GFF3.Record,1}, aln_interval::Interval{Bool}, args::Dict)
    """Process a single alignment interval for all valid overlaps with features."""
    aln_feat_overlaps = filter_alignment_feature_overlaps(features, aln_interval, isstranded(args["stranded"]))
    for feature in aln_feat_overlaps
        feature_name = get_feature_name_from_attrs(feature, args["attribute_type"])
        align_feat_ratio::Float32 = compute_align_feat_ratio(feat_overlaps[feature_name].coords_set, aln_interval)
        if align_feat_ratio > 0.0
            aln_type = gettype_alignment(metadata(aln_interval))
            increment_feature_overlap_information!(feat_overlaps[feature_name], align_feat_ratio, aln_type)
        end
    end
end

function process_overlaps!(feat_overlap::FeatureOverlap, multimapped_dict::Dict{String, IntervalCollection}, reader::BAM.Reader, feature::GFF3.Record, args::Dict)
    """Process all alignment intervals that overlap with feature intervals."""
    max_frag_size = convert(UInt, args["max_fragment_size"])
    gff_strand = getstrand(feature, isstranded(args["stranded"]))
    for record in eachoverlap(reader, feature, args["stranded"], max_frag_size)
        record_name = BAM.tempname(record)

        # Getting the interval here feels redundant since it is calculated in the Base.iterate function
        # But returning the record instead of the interval seemed to be much faster through initial testing
        aln_interval = get_alignment_interval(record, max_frag_size, is_reversestranded(args["stranded"]))
        @assert aln_interval !== nothing "All iterated records should be a validated read or fragment"
        aln_type = gettype_alignment(metadata(aln_interval))
        args["pp_only"] && isa(aln_type, ReadAlignment) && continue
        # Save multimapped records as Interval objects for later, if needed
        if ismultimapped(record)
            if !args["rm_multimap"]
                get!(multimapped_dict, record_name, IntervalCollection{Bool}())
                push!(multimapped_dict[record_name], aln_interval)
            end
            continue
        end
        process_overlap!(feat_overlap, aln_interval, gff_strand, args["stranded"])
    end
end

function process_overlap!(feat_overlap::FeatureOverlap, aln_interval::Interval{Bool}, gff_strand::Char, strand_type::String)
    """Process current single feature-alignment overlap."""
    bam_strand = getstrand(aln_interval, isstranded(strand_type))
    bam_strand == gff_strand || return
    align_feat_ratio::Float32 = compute_align_feat_ratio(feat_overlap.coords_set, aln_interval)
    if align_feat_ratio > 0.0
        aln_type = gettype_alignment(metadata(aln_interval))
        increment_feature_overlap_information!(feat_overlap, align_feat_ratio, aln_type)
    end
end
