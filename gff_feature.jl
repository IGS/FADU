#!/usr/bin/env julia

"""
gff_feature.jl - Houses functions that process GFF features.

By: Shaun Adkins (sadkins@som.umaryland.edu)
"""

isstranded(strand_type::String) = return (strand_type == "no" ? false : true)


function create_uniq_coords_dict(features::Array{GFF3.Record,1}, stranded_type::String)
    uniq_coords = Dict{String, Dict}()
    seqids = Set(GFF3.seqid(x) for x in features)
    for seqid in seqids
        uniq_coords[seqid] = get_nonoverlapping_coords_by_seqid(features, seqid)
        # If alignment data is unstranded, then symdiff both strands to + strand
        if !isstranded(stranded_type)
            # Ensure a default empty set is used if that strand has no coordinates
            symdiff!(get!(uniq_coords[seqid],'+', BitSet()), pop!(uniq_coords[seqid], '-', BitSet()))
        end
    end
    return uniq_coords
end

function get_feature_coords_set(feature::GFF3.Record)
    """Get the range of coordinates for this feature, returned as a Set."""
    return BitSet(GFF3.seqstart(feature) : GFF3.seqend(feature))
end

function get_featurename_from_attrs(feature::GFF3.Record, attr_type::String)
    """Get attribute ID to use as the feature name."""
    validate_feature_attribute(feature, attr_type)
    return GFF3.attributes(feature, attr_type)[1]
end

function get_feature_nonoverlapping_coords(feature::GFF3.Record, uniq_coords::Dict{String, Dict}, stranded::Bool)
    """Get set of nonoverlapping coordinates for given feature."""
    seqid = GFF3.seqid(feature)
    strand = getstrand(feature, stranded)
    feat_coords = get_feature_coords_set(feature)
    return intersect(feat_coords, uniq_coords[seqid][strand])
end

function get_nonoverlapping_coords_by_seqid(features::Array{GFF3.Record,1}, seqid::String)
    """ Get dictionary nonoverlapping coordinates based on the features on this sequence ID."""
    uniq_seqid_coords = Dict{Char,BitSet}()
    # Get only features with this reference id
    seqid_feats = filter(x -> GFF3.seqid(x) == seqid, features)
    strand = Set(map(x -> GFF3.strand(x), seqid_feats))
    for s in strand
        strandchar = convert(Char, s)
        get!(uniq_seqid_coords, strandchar, get_nonoverlapping_coords_by_seqid_and_strand(seqid_feats, s))
    end
    return uniq_seqid_coords
end

function get_nonoverlapping_coords_by_seqid_and_strand(features::Array{GFF3.Record,1}, strand::Strand)
    """ Get set of nonoverlapping coordinates based on the features for this strand of this sequence ID."""
    # Get features one strand at a time
    return BitSet(mapreduce((x -> GFF3.strand(x) == strand ? get_feature_coords_set(x) : BitSet()), symdiff, features))
end

function getstrand(feature::GFF3.Record, stranded::Bool)
    """Get strand of GFF3 feature with respect to strandedness arguments."""
    return getstrand(convert(Interval, feature), stranded)
end

function validate_feature_attribute(feature::GFF3.Record, attr_type::String)
    """Ensure only one attribute of this type exists for this feature record."""
    gene_vector = GFF3.attributes(feature, attr_type)    # col 9
    length(gene_vector) > 0 || error("ERROR - Attribute field 'ID' found to have no entries.\n")
    length(gene_vector) == 1 || error("ERROR - Attribute field 'ID' found to have multiple entries.\n")
end
