#!/usr/bin/env julia

"""
feature_counts.jl - Houses functions that deal with GFF feature-BAM alignment overlaps and counts

By: Shaun Adkins (sadkins@som.umaryland.edu)
"""


### Feature Overlap - Miscellaneous information about a feature with respect to its overlaps to alignments

mutable struct FeatureOverlap
    num_alignments::Float32
    feat_counts::Float32
    coords_set::BitSet
end

coordinate_set(f::FeatureOverlap) = f.coords_set
featurecounts(f::FeatureOverlap) = f.feat_counts
totalalignments(f::FeatureOverlap) = f.num_alignments

inc_featurecounts!(f::FeatureOverlap, counts) = f.feat_counts += counts
inc_alignments!(f::FeatureOverlap, num_alignments) = f.num_alignments += num_alignments

### MultimappedAlignment - Alignments that overlap features but are multimapped

struct MultimappedAlignment
    featurename::String
    align_feat_ratio::Float32
    alignmenttype::T where {T<:AbstractAlignment}
end

featurename(m::MultimappedAlignment) = m.featurename
alignment_feature_ratio(m::MultimappedAlignment) = m.align_feat_ratio
alignment_type(m::MultimappedAlignment) = m.alignmenttype


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
        totalcounts += featurecounts(feat_overlaps[featurename]) * count
    end
    return totalcounts
end

function calc_totalcounts(feat_overlaps::Dict{String, FeatureOverlap})
    """Calculate the sum of all the feature counts."""
    return sum(featurecounts(feat_overlaps[featurename]) for featurename in keys(feat_overlaps))
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

function compute_mm_counts_by_em(feat_overlaps::Dict{String, FeatureOverlap}, multimapped_dict::Dict{String,StructArray})
    """Adjust the feature counts of the multimapped overlaps via the Expectation-Maximization algorithm."""
    # Initialize multimapped overlap feature dictionary
    featurenames = collect(keys(feat_overlaps))
    adjusted_mm_overlaps = Dict{String, FeatureOverlap}(featurename => initialize_overlap_info(coordinate_set(feat_overlaps[featurename])) for featurename in featurenames)

    # Iterate through the multimaped alignments for each template
    @inbounds for (tn, templatealignments) in multimapped_dict
        # All multimapped alignments for a record should have the same alignment type
        alignmenttype = alignment_type(first(templatealignments))

        # I do not use "unique" for this because there may be a chance that the template aligns to the same feature multiple times
        overlappingfeatures = [featurename(a) for a in templatealignments]

        template_feature_counter = Dict{String, UInt}(fn => zero(UInt) for fn in unique(overlappingfeatures))
        template_feat_counts = Dict{String, Float32}(fn => zero(Float32) for fn in unique(overlappingfeatures))
        Threads.@threads for fn in unique(overlappingfeatures)
            # Filter all alignments that have this particular template-feature overlap
            templatefeaturealignments = filter(x -> x.featurename == fn, templatealignments)
            template_feature_counter[fn] += length(templatefeaturealignments)
            template_feat_counts[fn] += sum(alignment_feature_ratio, templatefeaturealignments)
        end

        template_totalcounts = calc_subset_totalcounts(feat_overlaps, template_feature_counter)
        # In cases where every feature for this template had no overlaps with any singly-mapped alignments, do not consider them.
        template_totalcounts > zero(Float32) || continue

        # Adjust contribution proportion of multimapped reads, based on estimated relative abundance
        for (featurename, counts) in template_feat_counts
            counts > zero(Float32) || continue
            adjusted_feat_counts = compute_mm_adjusted_counts(featurecounts(feat_overlaps[featurename]), counts, template_totalcounts)
            # If the feature had no contribution to any singly-mapped alignments, ignore it.
            adjusted_feat_counts > zero(Float32) && increment_feature_overlap_information!(adjusted_mm_overlaps[featurename], adjusted_feat_counts, alignmenttype)
        end
    end
    return adjusted_mm_overlaps
end

function create_feat_overlaps_dict(features::Array{GFF3.Record,1}, uniq_coords::Dict{String,Dict}, attr_type::String, stranded_type::String)
    feat_overlaps = Dict{String, FeatureOverlap}()
    for feature in features
        featurename = get_featurename_from_attrs(feature, attr_type)
        uniq_feat_coords = get_feature_nonoverlapping_coords(feature, uniq_coords, isstranded(stranded_type))
        feat_overlaps[featurename] = initialize_overlap_info(uniq_feat_coords)
    end
    return feat_overlaps
end

function filter_features_overlapping_alignments(features::Array{GFF3.Record,1}, alignment::Interval, isstranded::Bool)
    """Filter features to those just that align with the current alignment on the same strand."""
    alignmentstrand = getstrand(alignment, isstranded)
    overlapping_features = filter(x -> isoverlapping(convert(Interval, x), alignment), features)
    filter!(x -> alignmentstrand == getstrand(x, isstranded), overlapping_features)
    return overlapping_features   # Returns an array of GFF3 Records
end

function increment_feature_overlap_information!(feat_overlap::FeatureOverlap, align_feat_ratio::Float32, alignmenttype::T) where {T<:AbstractAlignment}
    """Increment number of alignments and feature count information for feature if alignment overlapped with uniq coords."""
    inc_alignments!(feat_overlap, multiplier(alignmenttype))
    inc_featurecounts!(feat_overlap, align_feat_ratio * multiplier(alignmenttype))
end

function initialize_overlap_info(uniq_feat_coords::BitSet)
    """Initialize overlap diction information."""
    num_alignments = zero(Float32); feat_counts = zero(Float32)
    return FeatureOverlap(num_alignments, feat_counts, uniq_feat_coords)
end

function merge_mm_counts(feat_overlaps::Dict{String, FeatureOverlap}, mm_feat_overlaps::Dict{String, FeatureOverlap}, adjust_num_alignments::Bool=true)
    """Merge EM-computed counts for multimapped reads and singly-mapped feature counts.  Return a new Dict of FeatureOverlaps"""
    adjusted_overlaps = deepcopy(feat_overlaps)
    featurenames = collect(keys(feat_overlaps))
    [inc_featurecounts!(adjusted_overlaps[fn], featurecounts(mm_feat_overlaps[fn])) for fn in featurenames]
    adjust_num_alignments || return adjusted_overlaps

    [inc_alignments!(adjusted_overlaps[fn], totalalignments(mm_feat_overlaps[fn])) for fn in featurenames]
    return adjusted_overlaps
end

function merge_mm_counts!(feat_overlaps::Dict{String, FeatureOverlap}, mm_feat_overlaps::Dict{String, FeatureOverlap}, adjust_num_alignments::Bool=true)
    """Merge EM-computed counts for multimapped reads back into the general feature overlap counts."""
    featurenames = collect(keys(feat_overlaps))
    [inc_featurecounts!(feat_overlaps[fn], featurecounts(mm_feat_overlaps[fn])) for fn in featurenames]
    adjust_num_alignments || return

    [inc_alignments!(feat_overlaps[fn], totalalignments(mm_feat_overlaps[fn])) for fn in featurenames]
end

function process_feature_overlaps!(feat_overlaps::Dict{String, FeatureOverlap}, multimapped_dict::Dict{String,StructArray}, reader::BAM.Reader, feature::GFF3.Record, args)
    """Process all alignment intervals that overlap with the given feature interval."""
    featurename = get_featurename_from_attrs(feature, args["attribute_type"])
    feat_overlap = feat_overlaps[featurename]
    max_frag_size = convert(UInt, args["max_fragment_size"])
    featurestrand = getstrand(feature, isstranded(args["stranded"]))
    pp_only = args["pp_only"]
    rm_multimap = args["rm_multimap"]

    # Returning SuperBAMRecord objects instead of BAM.Record objects
    for rec in eachoverlap(reader, feature, args["stranded"], max_frag_size)
        alignmentinterval = interval(rec)
        @assert alignmentinterval !== nothing "All iterated records should be a validated read or fragment"
        strand(rec) == featurestrand || continue

        alignmenttype = alignment_type(rec)
        pp_only && isa(alignmenttype, ReadAlignment) && continue

        align_feat_ratio =  compute_alignment_feature_ratio(coordinate_set(feat_overlap), alignmentinterval)::Float32
        # Save multimapped records as Interval objects for later, if needed, otherwise increment information for this feature now
        if ismultimapped(record(rec))
            if !rm_multimap
                templatename = BAM.tempname(record(rec))
                get!(multimapped_dict, templatename, StructArray{MultimappedAlignment}(undef,0))
                push!(multimapped_dict[templatename], MultimappedAlignment(featurename, align_feat_ratio, alignmenttype))
            end
            continue
        end
        align_feat_ratio > 0.0 && increment_feature_overlap_information!(feat_overlap, align_feat_ratio, alignmenttype)
    end

    # NOTE: This function modifies feat_overlaps and multimapped_dict
end

function process_all_feature_overlaps(feat_overlaps::Dict{String,FeatureOverlap}, features::Vector{GFF3.Record}, reader::BAM.Reader, args)
    """Process all features from the annotation file for overlaps with BAM alignments."""
    multimapped_dict = Dict{String,StructArray}()
    [process_feature_overlaps!(feat_overlaps, multimapped_dict, reader, feature, args) for feature in features]
    return multimapped_dict
end
