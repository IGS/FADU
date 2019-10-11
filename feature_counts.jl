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

function calc_relative_abundance(feature_counts::Float32, template_totalcounts::Float32)::Float32
    """Calculate relative abundance fraction for the given feature with respect to all features overlapping the template."""
    @assert template_totalcounts > 0 "The template_totalcounts value should be greater than 0"
    # Total relative abundances for all features overlapping this template should = 1.0
    return @fastmath feature_counts / template_totalcounts
end

function calc_subset_totalcounts(feat_overlaps::Dict{String, FeatureOverlap}, subset_feat_overlaps::Array{GFF3.Record,1}, attribute_type::String)
    """For all features overlapping the alignment template, return the sum of all those feature counts."""
    totalcounts = zero(Float32)
    # Worth noting if a feature overlaps multiple times, it will be reflected in totalcounts
    # This will be appropriately handled when all alignments are iterated through to increment the FeatureOverlap feat_counts for the feature.
    for feature in subset_feat_overlaps
        featurename = get_featurename_from_attrs(feature, attribute_type)
        totalcounts += feat_overlaps[featurename].feat_counts
    end
    return totalcounts
end

function calc_totalcounts(feat_overlaps::Dict{String, FeatureOverlap})
    """Calculate the sum of all the feature counts."""
    return sum(feat_overlaps[featurename].feat_counts for featurename in keys(feat_overlaps))
end

function calc_tpm(len::UInt, totalcounts::Float32, feat_counts::Float32)
    """Calculate TPM score for current feature."""
    return @fastmath(feat_counts * 1000 / len) * 1000000 / totalcounts
end

function compute_alignment_feature_ratio(uniq_feat_coords::Set{UInt}, alignment::Interval)::Float32
    """Calculate the ratio of fragament coordinates that intersect with non-overlapping feature coordinates."""
    # Pertinent alignment info
    alignment_coords =  get_alignment_coords_set(alignment)
    alignment_intersect = intersect(uniq_feat_coords, alignment_coords)
    # Percentage of alignment that aligned with this annotation feature
    # Return as Float32 for performance
    return @fastmath length(alignment_intersect) / length(alignment_coords)
end

function compute_mm_counts_by_em(uniq_feat_overlaps::Dict{String, FeatureOverlap}, mm_feat_overlaps::Dict{String, FeatureOverlap}, alignment_dict::Dict{String,IntervalCollection}, features::Array{GFF3.Record,1}, args::Dict)
    """Adjust the feature counts of the multimapped overlaps via the Expectation-Maximization algorithm."""
    #TODO: Clean up function
    adjusted_mm_overlaps = Dict{String, FeatureOverlap}()
    uniq_and_mm_overlaps = deepcopy(uniq_feat_overlaps)
    for featurename in keys(uniq_feat_overlaps)
        adjusted_mm_overlaps[featurename] = initialize_overlap_info(uniq_feat_overlaps[featurename].coords_set)
        # Combine singly-mapped and multimapped feature counts
        # Reason is to eventually update relative abundances per feature per alignment overlap
        uniq_and_mm_overlaps[featurename].feat_counts += mm_feat_overlaps[featurename].feat_counts
    end

    for record_tempname in keys(alignment_dict)
        # Get all features that align with the given multimapped template name on the same strand
        template_features = Array{GenomicFeatures.GFF3.Record,1}()
        for aln_interval in alignment_dict[record_tempname]
            append!(template_features, filter_alignment_feature_overlaps(features, aln_interval, isstranded(args["stranded"])))
        end

        template_totalcounts = calc_subset_totalcounts(uniq_and_mm_overlaps, template_features, args["attribute_type"])
        # In cases where every feature for this template had no overlaps with any singly-mapped alignments, do not consider them.
        template_totalcounts > 0 || continue

        # Adjust contribution proportion of multimapped reads, based on estimated relative abundance
        for aln_interval in alignment_dict[record_tempname]
            aln_feat_overlaps = filter_alignment_feature_overlaps(features, aln_interval, isstranded(args["stranded"]))
            for feature_overlap in aln_feat_overlaps
                featurename = get_featurename_from_attrs(feature_overlap, args["attribute_type"])
                featurestrand = getstrand(feature_overlap, isstranded(args["stranded"]))
                relative_abundance = calc_relative_abundance(uniq_and_mm_overlaps[featurename].feat_counts, template_totalcounts)
                adjusted_overlap = adjusted_mm_overlaps[featurename]
                process_overlap!(adjusted_overlap, aln_interval, relative_abundance)
            end
        end
    end
    return adjusted_mm_overlaps
end

function create_feat_overlaps_dict(features::Array{GenomicFeatures.GFF3.Record,1}, uniq_coords::Dict{String,Dict}, attr_type::String, stranded_type::String)
    feat_overlaps = Dict{String, FeatureOverlap}()
    for feature in features
        featurename = get_featurename_from_attrs(feature, attr_type)
        uniq_feat_coords = get_feature_nonoverlapping_coords(feature, uniq_coords, isstranded(stranded_type))
        feat_overlaps[featurename] = initialize_overlap_info(uniq_feat_coords)
    end
    return feat_overlaps
end

function filter_alignment_feature_overlaps(features::Array{GenomicFeatures.GFF3.Record,1}, alignment::Interval, isstranded::Bool)
    """Filter features to those just that align with the current alignment on the same strand."""
    aln_strand = getstrand(alignment, isstranded)
    aln_feat_overlaps = filter(x -> isoverlapping(convert(Interval, x), alignment), features)
    filter!(x -> aln_strand == getstrand(x, isstranded), aln_feat_overlaps)
    return aln_feat_overlaps
end

function increment_feature_overlap_information!(feat_overlap::FeatureOverlap, align_feat_ratio::Float32, aln_type::T) where {T<:AbstractAlignment}
    """Increment number of alignments and feature count information for feature if alignment overlapped with uniq coords."""
    feat_overlap.num_alignments += aln_type.count_multiplier
    feat_overlap.feat_counts += align_feat_ratio * aln_type.count_multiplier
end

function initialize_overlap_info(uniq_feat_coords::Set{UInt})
    """Initialize overlap diction information."""
    num_alignments = zero(Float32); feat_counts = zero(Float32)
    return FeatureOverlap(num_alignments, feat_counts, uniq_feat_coords)
end

function merge_mm_counts!(feat_overlaps::Dict{String, FeatureOverlap}, mm_feat_overlaps::Dict{String, FeatureOverlap})
    """Merge EM-computed counts for multimapped reads back into the general feature overlap counts."""
    for featurename in keys(feat_overlaps)
        feat_overlaps[featurename].num_alignments += mm_feat_overlaps[featurename].num_alignments
        feat_overlaps[featurename].feat_counts += mm_feat_overlaps[featurename].feat_counts
    end
end

function process_overlaps!(feat_overlap::FeatureOverlap, multimapped_dict::Dict{String, IntervalCollection}, reader::BAM.Reader, feature::GFF3.Record, args::Dict)
    """Process all alignment intervals that overlap with feature intervals."""
    max_frag_size = convert(UInt, args["max_fragment_size"])
    featurestrand = getstrand(feature, isstranded(args["stranded"]))
    for record in eachoverlap(reader, feature, args["stranded"], max_frag_size)
        # Getting the interval here feels redundant since it is calculated in the Base.iterate function
        # But returning the record instead of the interval seemed to be much faster through initial testing
        aln_interval = get_alignment_interval(record, max_frag_size, is_reversestranded(args["stranded"]))
        alignmentstrand = getstrand(aln_interval, isstranded(args["stranded"]))
        alignmentstrand == featurestrand || continue
        @assert aln_interval !== nothing "All iterated records should be a validated read or fragment"
        aln_type = gettype_alignment(metadata(aln_interval))
        args["pp_only"] && isa(aln_type, ReadAlignment) && continue
        # Save multimapped records as Interval objects for later, if needed
        if ismultimapped(record)
            if !args["rm_multimap"]
                template_name = BAM.tempname(record)
                get!(multimapped_dict, template_name, IntervalCollection{Bool}())
                push!(multimapped_dict[template_name], aln_interval)
            end
            continue
        end
        process_overlap!(feat_overlap, aln_interval)
    end
end

function process_overlap!(feat_overlap::FeatureOverlap, aln_interval::Interval{Bool}, relative_abundance::Float32=1.0f0)
    """Process current single feature-alignment overlap."""
    align_feat_ratio::Float32 = compute_alignment_feature_ratio(feat_overlap.coords_set, aln_interval)
    if align_feat_ratio > 0.0
        aln_type = gettype_alignment(metadata(aln_interval))
        adjusted_align_feat_ratio = align_feat_ratio * relative_abundance
        increment_feature_overlap_information!(feat_overlap, adjusted_align_feat_ratio, aln_type)
    end
end
