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

function calc_tpm(len::UInt, depth_sum::Float32, feat_counts::Float32)
    """Calculate TPM score for current feature."""
    return @fastmath(feat_counts *1000 / len) * 1000000 / depth_sum
end

function calc_total_counts(feat_overlaps::Dict{String, FeatureOverlap})
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