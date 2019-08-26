#!/usr/bin/env julia

"""
bam_record.jl - Houses functions that deal with BAM Record objects and properties

By: Shaun Adkins (sadkins@som.umaryland.edu)
"""

abstract type AbstractAlignment end

struct FragmentAlignment<:AbstractAlignment 
    count_multiplier::Float32
end
struct ReadAlignment<:AbstractAlignment 
    count_multiplier::Float32
end

FragmentAlignment() = FragmentAlignment(1.0)
ReadAlignment() = ReadAlignment(0.5)

is_reverse_stranded(strand_type::String) = return (strand_type == "reverse" ? true : false)

#isduplicate(record::BAM.Record) = BAM.flag(record) & SAM.FLAG_DUP == 0x0400
#ismatereverse(record::BAM.Record) = BAM.flag(record) & SAM.FLAG_MREVERSE == 0x0020
isproperpair(record::BAM.Record) = BAM.flag(record) & SAM.FLAG_PROPER_PAIR == 0x0002
isread1(record::BAM.Record) = BAM.flag(record) & SAM.FLAG_READ1 == 0x0040
isread2(record::BAM.Record) = BAM.flag(record) & SAM.FLAG_READ2 == 0x0080
isreverse(record::BAM.Record) = BAM.flag(record) & SAM.FLAG_REVERSE == 0x0010

# Forward-stranded assay
## Positive Strand:
### R1 - forward (65), R2 - revcom (145)
## Negative Strand
### R1 - revcom (81), R2 - forward (129)
# Reverse-stranded assay (i.e. Illumina) have these strands flipped

### For non-properly paired reads
# Singletons will either have flag 0 or 16

# Forward strand flags
flag65(record::BAM.Record) = isread1(record) && !isreverse(record)
flag145(record::BAM.Record) = isread2(record) && isreverse(record)
# Reverse strand flags
flag81(record::BAM.Record) = isread1(record) && isreverse(record)
flag129(record::BAM.Record) = isread2(record) && !isreverse(record)

function assign_read_to_strand(record::BAM.Record, reverse_strand::Bool=false)
    """Use the bitwise flags to assign the paired read to the correct strand."""

    pos_flags = Set{Bool}([flag145(record), flag65(record)])
    neg_flags = Set{Bool}([flag81(record), flag129(record)])
    if reverse_strand
        pos_flags = Set{Bool}([flag81(record), flag129(record)])
        neg_flags = Set{Bool}([flag145(record), flag65(record)])
    end

    # Does record flag belong in the pos or neg set?
    any(pos_flags) && return '+'
    any(neg_flags) && return '-'

    # Read must be a singleton
    if isreverse(record)
        # If reverse-stranded, strand is flipped
        reverse_strand && return '+'
        return '-'
    end
    reverse_strand && return '-'
    return '+'
end

#function determine_record_type(record::BAM.Record, max_frag_size::UInt)
function canbefragment(record::BAM.Record, max_frag_size::UInt)
    """Determine if BAM Record is a read (R) or fragment (F)."""
    # True = Can validate as a fragment
    # False = Can validate as a read
    if validate_fragment(record) && is_templength_smaller_than_max_fragment_size(BAM.templength(record), max_frag_size)
        return true
        #return FragmentAlignment()
    elseif validate_read(record, max_frag_size)
        return false
        #return ReadAlignment()
    end
    # For fragments, only one read is looked at.  The other one is essentially skipped to avoid overcounting.
    # Reads that also fail validation go here.
    return nothing
end

function get_alignment_start_end(record::BAM.Record, isfragment::Bool)
    return get_alignment_start_end(record, gettype_alignment(isfragment))
end

function get_alignment_start_end(record::BAM.Record, record_type::FragmentAlignment)
    """Get the start and end coordinates of the fragment.  Assuming the use of the negative template"""
    return BAM.nextposition(record):BAM.rightposition(record)
end

function get_alignment_start_end(record::BAM.Record, record_type::ReadAlignment)
    """Get the start and end coordinates of the read."""
    return BAM.position(record):BAM.rightposition(record)
end

function get_alignment_interval(record::BAM.Record, max_frag_size::UInt, reverse_strand::Bool=false)
    """Return an alignment-based Interval for the current record.""" 
    #record_type = determine_record_type(record, max_frag_size)
    #record_type === nothing && return nothing  # Skip alignments not being considered
    isfragment = canbefragment(record, max_frag_size)
    isfragment === nothing && return nothing  # Skip alignments not being considered
    return Interval(BAM.refname(record),
                    get_alignment_start_end(record, isfragment),
                    assign_read_to_strand(record, reverse_strand), 
                    isfragment
    )
end

function get_alignment_coords_set(alignment::Interval)
    """Get the range of coordinates for this alignment, returned as a Set."""
    return Set{UInt}(leftposition(alignment) : rightposition(alignment))
end

function get_alignment_interval(record::BAM.Record, max_frag_size::UInt, strand_type::String)
    return get_alignment_interval(record, max_frag_size, is_reverse_stranded(strand_type))
end

function getstrand(interval::Interval, stranded::Bool)
    """Get strand of interval with respect to strandedness arguments."""
    strand = '+'
    if stranded
        strand = convert(Char, interval.strand)
    end
    return strand
end

function gettype_alignment(isfragment::Bool)
    return (isfragment ? FragmentAlignment() : ReadAlignment())
end

function ismultimapped(record::BAM.Record)
    """Test to see if alignment is multimapped across multiple regions of the genome."""
    try
        return haskey(record, "NH") && record["NH"] > 1
    catch
        # Ran into bug where some attributes are 'X' type which is not valid
        # Bioalignments BAM.Record.auxdata(record) throws LoadError for these
        return false
    end
end

function is_templength_negative(templength::Int64)
    """Check to see if the read template is going in the opposite version."""
    return templength < 0
end

function is_templength_smaller_than_max_fragment_size(templength::Int64, max_frag_size::UInt)
    """Check to see if fragment template length is smaller than specified maximum fragment size."""
    return abs(templength) <= max_frag_size
end

function validate_fragment(record::BAM.Record)
    """Ensure alignment record can be used to calculate fragment depth. (Only need negative template of fragment)"""
    return isproperpair(record) && is_templength_negative(BAM.templength(record))
end

function validate_read(record::BAM.Record, max_frag_size::UInt)
    """Ensure alignment record can be used to calculate read depth."""
    return !(isproperpair(record) && is_templength_smaller_than_max_fragment_size(BAM.templength(record), max_frag_size))
end
