#!/usr/bin/env julia

"""
gff_feature.jl - Houses functions that process GFF features.

By: Shaun Adkins (sadkins@som.umaryland.edu)
"""

isstranded(strand_type::String) = return (strand_type == "no" ? false : true)


function create_uniq_coords_dict(features::Array{GenomicFeatures.GFF3.Record,1}, stranded_type::String)
    uniq_coords = Dict{String, Dict}()
    seqids = Set(map(x -> GFF3.seqid(x), features))
    for seqid in seqids
        uniq_coords[seqid] = get_nonoverlapping_coords_by_seqid(features, seqid)
        # If alignment data is unstranded, then symdiff both strands to + strand 
        if !isstranded(stranded_type)
            symdiff!(uniq_coords[seqid]['+'], pop!(uniq_coords[seqid], '-', Set()))
        end
    end
    return uniq_coords
end

function get_feature_coords_set(feature::GenomicFeatures.GFF3.Record)
    """Get the range of coordinates for this feature, returned as a Set."""
    return Set{UInt}(leftposition(feature) : rightposition(feature))
end

function get_feature_name_from_attrs(feature::GFF3.Record, attr_type::String)
    """Get attribute ID to use as the feature name."""
    gene_vector = GFF3.attributes(feature, attr_type)    # col 9
    validate_feature_attribute(gene_vector)
    return gene_vector[1]
end

function get_feature_nonoverlapping_coords(feature::GenomicFeatures.GFF3.Record, uniq_coords::Dict{String, Dict}, stranded::Bool)
    """Get set of nonoverlapping coordinates for given feature."""
    seqid = GFF3.seqid(feature)
    strand = getstrand(feature, stranded)
    feat_coords = get_feature_coords_set(feature)
    return intersect(feat_coords, uniq_coords[seqid][strand])
end

function get_nonoverlapping_coords_by_seqid(features::Array{GenomicFeatures.GFF3.Record,1}, seqid::String)
    """ Get dictionary nonoverlapping coordinates based on the features on this sequence ID."""
    uniq_seqid_coords = Dict{Char,Set{UInt}}()
    # Get only features with this reference id
    seqid_feats = filter(x -> GFF3.seqid(x) == seqid, features)
    strand = Set(map(x -> GFF3.strand(x), seqid_feats))
    for s in strand
        strandchar = convert(Char, s)
        uniq_seqid_coords[strandchar] = get_nonoverlapping_coords_by_seqid_and_strand(seqid_feats, s)
    end
    return uniq_seqid_coords
end

function get_nonoverlapping_coords_by_seqid_and_strand(features::Array{GenomicFeatures.GFF3.Record,1}, strand::Strand)
    """ Get set of nonoverlapping coordinates based on the features for this strand of this sequence ID."""
    # Get features one strand at a time
    seqid_feats_by_strand = filter(x-> GFF3.strand(x) == strand, features)
    return Set{UInt}(mapreduce(x -> get_feature_coords_set(x), symdiff, seqid_feats_by_strand))
end

function getstrand(feature::GenomicFeatures.GFF3.Record, stranded::Bool)
    """Get strand of GFF3 feature with respect to strandedness arguments."""
    return getstrand(convert(Interval, feature), stranded)
end

function validate_feature_attribute(gene_vector::Vector{String})
    """Ensure only one attribute of this type exists for this feature record."""
    length(gene_vector) > 0 || error("ERROR - Attribute field 'ID' found to have no entries.\n")
    length(gene_vector) == 1 || error("ERROR - Attribute field 'ID' found to have multiple entries.\n")
end