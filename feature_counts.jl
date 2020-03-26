#!/usr/bin/env julia

"""
feature_counts.jl - Houses functions that deal with GFF feature-BAM alignment overlaps and counts

By: Shaun Adkins (sadkins@som.umaryland.edu)
"""

mutable struct FeatureOverlap
    num_alignments::Float32
    feat_counts::Float32
    coords_set::BitSet
end

function calc_relative_abundance(feature_counts::Float32, template_totalcounts::Float32)::Float32
    """Calculate relative abundance fraction for the given feature with respect to all features overlapping the template."""
    @assert template_totalcounts > 0 "The template_totalcounts value should be greater than 0"
    # Total relative abundances for all features overlapping this template should = 1.0
    return @fastmath feature_counts / template_totalcounts
end

function calc_subset_totalcounts(feat_overlaps::Dict{String, FeatureOverlap}, feature_counter::Dict{String, UInt})
    """For all features overlapping the alignment template, return the sum of all those feature counts."""
    totalcounts = zero(Float32)
    # Worth noting if a feature overlaps multiple times, it will be reflected in totalcounts
    # This will be appropriately handled when all alignments are iterated through to increment the FeatureOverlap feat_counts for the feature.
    for (featurename, count) in feature_counter
        totalcounts += feat_overlaps[featurename].feat_counts * count
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

function compute_alignment_feature_ratio(uniq_feat_coords::BitSet, alignment::Interval)::Float32
    """Calculate the ratio of fragament coordinates that intersect with non-overlapping feature coordinates."""
    # Pertinent alignment info
    alignment_coords =  get_alignment_coords_set(alignment)
    alignment_intersect = intersect(uniq_feat_coords, alignment_coords)
    # Percentage of alignment that aligned with this annotation feature
    # Return as Float32 for performance
    return @fastmath length(alignment_intersect) / length(alignment_coords)
end

function compute_mm_adjusted_counts(total_feat_counts::Float32, template_feat_counts::Float32, template_totalcounts::Float32)::Float32
    """Calculate the adjustment feature counts for a multimapped template to a given feature."""
    relative_abundance = calc_relative_abundance(total_feat_counts, template_totalcounts)
    return template_feat_counts * relative_abundance
end

function compute_mm_counts_by_em(feat_overlaps::Dict{String, FeatureOverlap}, alignment_dict::Dict{String,IntervalCollection}, features::Array{GFF3.Record,1}, args::Dict)
    """Adjust the feature counts of the multimapped overlaps via the Expectation-Maximization algorithm."""
    # Initialize multimapped overlap feature dictionary
    featurenames = collect(keys(feat_overlaps))
    adjusted_mm_overlaps = Dict{String, FeatureOverlap}(featurename => initialize_overlap_info(feat_overlaps[featurename].coords_set) for featurename in featurenames)

    for record_tempname in keys(alignment_dict)
        # Get all features that align with the given multimapped template name on the same strand

        template_featurename_counter = Dict{String, UInt}(featurename => zero(UInt) for featurename in featurenames)
        template_feat_counts = Dict{String, Float32}(featurename => zero(Float32) for featurename in featurenames)

        # All multimapped alignments for a record should have the same alignment type
        alignment_type = getalignmenttype(metadata(first(alignment_dict[record_tempname])))

        for aln_interval in alignment_dict[record_tempname]
            overlapping_features = filter_features_overlapping_alignments(features, aln_interval, isstranded(args["stranded"]))
            # Only deal with alignments that overlap with a feature
            isempty(overlapping_features) && continue
            # Process each feature found in the alignment
            for feature in overlapping_features
                featurename = get_featurename_from_attrs(feature, args["attribute_type"])
                template_featurename_counter[featurename] += one(UInt)
                align_feat_ratio = compute_alignment_feature_ratio(feat_overlaps[featurename].coords_set, aln_interval)
                # Some features may not have contributions to any singly-mapped alignments.  They will not factor into the adjustments later
                template_feat_counts[featurename] += align_feat_ratio
            end
        end

        template_totalcounts = calc_subset_totalcounts(feat_overlaps, template_featurename_counter)
        # In cases where every feature for this template had no overlaps with any singly-mapped alignments, do not consider them.
        template_totalcounts > zero(Float32) || continue

        # Adjust contribution proportion of multimapped reads, based on estimated relative abundance
        for (featurename, counts) in template_feat_counts
            counts > zero(Float32) || continue
            adjusted_feat_counts = compute_mm_adjusted_counts(feat_overlaps[featurename].feat_counts, counts, template_totalcounts)
            # If the feature had no contribution to any singly-mapped alignments, ignore it.
            adjusted_feat_counts > zero(Float32) && increment_feature_overlap_information!(adjusted_mm_overlaps[featurename], adjusted_feat_counts, alignment_type)
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

function filter_features_overlapping_alignments(features::Array{GenomicFeatures.GFF3.Record,1}, alignment::Interval, isstranded::Bool)
    """Filter features to those just that align with the current alignment on the same strand."""
    alignmentstrand = getstrand(alignment, isstranded)
    overlapping_features = filter(x -> isoverlapping(convert(Interval, x), alignment), features)
    filter!(x -> alignmentstrand == getstrand(x, isstranded), overlapping_features)
    return overlapping_features   # Returns an array of GFF3 Records
end

function increment_feature_overlap_information!(feat_overlap::FeatureOverlap, align_feat_ratio::Float32, alignment_type::T) where {T<:AbstractAlignment}
    """Increment number of alignments and feature count information for feature if alignment overlapped with uniq coords."""
    feat_overlap.num_alignments += alignment_type.count_multiplier
    feat_overlap.feat_counts += align_feat_ratio * alignment_type.count_multiplier
end

function initialize_overlap_info(uniq_feat_coords::BitSet)
    """Initialize overlap diction information."""
    num_alignments = zero(Float32); feat_counts = zero(Float32)
    return FeatureOverlap(num_alignments, feat_counts, uniq_feat_coords)
end

function merge_mm_counts(feat_overlaps::Dict{String, FeatureOverlap}, mm_feat_overlaps::Dict{String, FeatureOverlap}, adjust_num_alignments::Bool=true)
    """Merge EM-computed counts for multimapped reads and singly-mapped feature counts.  Return a new Dict of FeatureOverlaps"""
    adjusted_overlaps = deepcopy(feat_overlaps)
    for featurename in keys(adjusted_overlaps)
        if adjust_num_alignments
            adjusted_overlaps[featurename].num_alignments += mm_feat_overlaps[featurename].num_alignments
        end
        adjusted_overlaps[featurename].feat_counts += mm_feat_overlaps[featurename].feat_counts
    end
    return adjusted_overlaps
end

function merge_mm_counts!(feat_overlaps::Dict{String, FeatureOverlap}, mm_feat_overlaps::Dict{String, FeatureOverlap}, adjust_num_alignments::Bool=true)
    """Merge EM-computed counts for multimapped reads back into the general feature overlap counts."""
    for featurename in keys(feat_overlaps)
        if adjust_num_alignments
            feat_overlaps[featurename].num_alignments += mm_feat_overlaps[featurename].num_alignments
        end
        feat_overlaps[featurename].feat_counts += mm_feat_overlaps[featurename].feat_counts
    end
end

# function process_overlaps!(feat_overlap::Dict{String, FeatureOverlap}, multimapped_dict::Dict{String, IntervalCollection}, reader::BAM.Reader, features::Array{GenomicFeatures.GFF3.Record,1}, args::Dict)
#     """Process all alignment intervals that overlap with feature intervals."""
#     record = BAM.Record()
#     max_frag_size = convert(UInt, args["max_fragment_size"])

#     while !eof(reader)
#         read!(reader, record)
#         # The following steps are figuring out if the record is worth looking at
#         # 1. Is the record mapped to the reference annotation?
#         BAM.ismapped(record) || continue

#         # 2. Establish alignment-based information.  Does the record pass validation as a read or a fragment we can use?
#         alignmentinterval = get_alignment_interval(record, max_frag_size, is_reversestranded(args["stranded"]))
#         alignmentinterval === nothing && continue

#         # 3. If only reading fragments (properly-paired reads), skip all non-properly-paired reads
#         alignmenttype = getalignmenttype(metadata(alignmentinterval))
#         args["pp_only"] && isa(alignmenttype, ReadAlignment) && continue

#         # 4. Only deal with alignments that overlap with a feature
#         overlapping_features = filter_features_overlapping_alignments(features, alignmentinterval, isstranded(args["stranded"]))
#         isempty(overlapping_features) && continue

#         # 5. Save multimapped records as Interval objects for later, if needed
#         if ismultimapped(record)
#             if !args["rm_multimap"]
#                 templatename = BAM.tempname(record)
#                 get!(multimapped_dict, templatename, IntervalCollection{Bool}())
#                 push!(multimapped_dict[templatename], alignmentinterval)
#             end
#             continue
#         end

#         # 6. Read/Fragment alignment can now be compared with overlapping part of feature.
#         for feature in overlapping_features
#             featurename = get_featurename_from_attrs(feature, args["attribute_type"])
#             process_overlap!(feat_overlap[featurename], alignmentinterval)
#         end
#     end
# end

# function process_overlap!(feat_overlap::FeatureOverlap, alignmentinterval::Interval{Bool})
#     """Process current single feature-alignment overlap."""
#     align_feat_ratio::Float32 = compute_alignment_feature_ratio(feat_overlap.coords_set, alignmentinterval)
#     if align_feat_ratio > zero(Float32)
#         alignment_type = getalignmenttype(metadata(alignmentinterval))
#         increment_feature_overlap_information!(feat_overlap, align_feat_ratio, alignment_type)
#     end
# end

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
        aln_type = getalignmenttype(metadata(aln_interval))
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

function process_overlap!(feat_overlap::FeatureOverlap, aln_interval::Interval{Bool})
    """Process current single feature-alignment overlap."""
    align_feat_ratio::Float32 = compute_alignment_feature_ratio(feat_overlap.coords_set, aln_interval)
    if align_feat_ratio > 0.0
        alignment_type = getalignmenttype(metadata(aln_interval))
        increment_feature_overlap_information!(feat_overlap, align_feat_ratio, alignment_type)
    end
end
