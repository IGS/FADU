#!/usr/bin/env julia

"""
alignment_overlaps.jl - An alternative implemention of the BioAlignments.jl
BAM overlap.jl code to convert Record types into fragment-based and read-based Interval types

Original code at https://github.com/BioJulia/BioAlignments.jl/blob/master/src/bam/overlap.jl

By: Shaun Adkins (sadkins@som.umaryland.edu)
"""

import BGZFStreams

### OverlapIterator - Extended version of original BAM OverlapIterator

struct OverlapIterator{T}
    reader::BAM.Reader{T}
    refname::String
    interval::UnitRange{Int}
    strand_type::String
    max_frag_size::UInt
end

function GenomicFeatures.eachoverlap(reader::BAM.Reader, interval::Interval, strand_type::String="no", max_frag_size::UInt=1000)
    return GenomicFeatures.eachoverlap(reader, interval.seqname, interval.first:interval.last, strand_type, max_frag_size)
end

function GenomicFeatures.eachoverlap(reader::BAM.Reader, interval, strand_type::String="no", max_frag_size::UInt=1000)
    return GenomicFeatures.eachoverlap(reader, convert(Interval, interval), strand_type, max_frag_size)
end

function GenomicFeatures.eachoverlap(reader::BAM.Reader, refname::AbstractString, interval::UnitRange, strand_type::String="no", max_frag_size::UInt=1000)
    return OverlapIterator(reader, String(refname), interval, strand_type, max_frag_size)
end

### SuperBAMRecord - BAM Record with additional information

struct SuperBAMRecord
    record::BAM.Record
    interval::Interval{Bool}
    strand::Char
    #alignmenttype::T where {T<:AbstractAlignment}
end

#function SuperBAMRecord(record::BAM.Record, interval::Interval, strand::Char)
#    return SuperBAMRecord(record, interval, strand, getalignmenttype(metadata(interval)))
#end

record(s::SuperBAMRecord) = s.record
interval(s::SuperBAMRecord) = s.interval
strand(s::SuperBAMRecord) = s.strand
alignment_type(s::SuperBAMRecord) = getalignmenttype(metadata(interval(s)))

# Iterator
# --------

function Base.iterate(iter::OverlapIterator)
    refindex = findfirst(isequal(iter.refname), iter.reader.refseqnames)
    if refindex == 0
        throw(ArgumentError("sequence name $(iter.refname) is not found in the header"))
    end
    @assert iter.reader.index !== nothing
    chunks = Indexes.overlapchunks(iter.reader.index.index, refindex, iter.interval)
    if !isempty(chunks)
        seek(iter.reader, first(chunks).start)
    end
    state = BAM.OverlapIteratorState(refindex, chunks, 1, BAM.Record())
    return iterate(iter, state)
end

function Base.iterate(iter::OverlapIterator, state)
    while state.chunkid ≤ lastindex(state.chunks)
        chunk = state.chunks[state.chunkid]

        # Ran into issue where the BAM.Reader processed all records but kept going which 
        # lead to the stream block_index to exceed the number of blocks.  This solves that.
        if eof(iter.reader.stream)
            break
        end

        while BGZFStreams.virtualoffset(iter.reader.stream) < chunk.stop
            read!(iter.reader, state.record)
            # Determine if fragment or read alignment.  Fragment alignments may overlap the interval whereas the single read may be to the right
            alignmentinterval = get_alignment_interval(state.record, iter.max_frag_size, iter.strand_type)
            alignmentinterval === nothing && continue

            c = GenomicFeatures.compare_overlap(alignmentinterval, Interval(iter.refname, iter.interval), isless)
            if c == 0
                alignmentstrand = getstrand(alignmentinterval, isstranded(iter.strand_type))
                return SuperBAMRecord(copy(state.record), alignmentinterval, alignmentstrand), state
            elseif c > 0
                # no more overlapping records in this chunk since records are sorted
                break
            end

        end
        state.chunkid += 1
        if state.chunkid ≤ lastindex(state.chunks)
            seek(iter.reader, state.chunks[state.chunkid].start)
        end
    end
    return nothing
end
