#!/usr/bin/env julia

"""
bam_record.jl - Houses functions that process BAM Record objects.

By: Shaun Adkins (sadkins@som.umaryland.edu)
"""

#is_duplicate(record::BAM.Record) = BAM.flag(record) & SAM.FLAG_DUP == 0x0400
#is_mate_reverse(record::BAM.Record) = BAM.flag(record) & SAM.FLAG_MREVERSE == 0x0020
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

function determine_record_type(record::BAM.Record, max_frag_size::UInt, pp_only::Bool)
    """Determine if BAM Record is a read (R) or fragment (F)."""
    if validate_fragment(record) && is_templength_smaller_than_max_fragment_size(BAM.templength(record), max_frag_size)
        return 'F'
    elseif !pp_only && validate_read(record, max_frag_size)
        return 'R'
    end
    # For fragments, only one read is looked at.  The other one is essentially skipped to avoid overcounting.
    # Reads that also fail validation go here.
    return nothing
end

function get_alignment_start_end(record::BAM.Record, record_type::Char)
    """Get the start and end coordinates of the fragment/read."""
    if record_type == 'R'
        return BAM.position(record):BAM.rightposition(record)
    end
    return BAM.nextposition(record):BAM.rightposition(record)
end

function get_alignment_interval(record::BAM.Record, record_type::Char, reverse_strand::Bool=false)
    """Return an alignment-based Interval for the current record.""" 
    return Interval(BAM.refname(record),
                    get_alignment_start_end(record, record_type),
                    assign_read_to_strand(record, reverse_strand), 
                    record_type
    )
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

function validate_fragment(record::BAM.Record)
    """Ensure alignment record can be used to calculate fragment depth. (Only need negative template of fragment)"""
    return is_proper_pair(record) && is_templength_negative(BAM.templength(record))
end

function validate_read(record::BAM.Record, max_frag_size::UInt)
    """Ensure alignment record can be used to calculate read depth."""
    return !(is_proper_pair(record) && is_templength_smaller_than_max_fragment_size(BAM.templength(record), max_frag_size))
end