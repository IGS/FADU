#!/usr/bin/env julia

"""
alignment_overlaps.jl - An alternative implemention of the BioAlignments.jl
BAM overlap.jl code to convert Record types into fragment-based and read-based Interval types

Original code at https://github.com/BioJulia/BioAlignments.jl/blob/master/src/bam/overlap.jl

By: Shaun Adkins (sadkins@som.umaryland.edu)
"""

import BGZFStreams

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

# Iterator
# --------

function Base.iterate(iter::OverlapIterator)
    refindex = findfirst(isequal(iter.refname), iter.reader.refseqnames)
    if refindex == 0
        throw(ArgumentError("sequence name $(iter.refname) is not found in the header"))
    end
    @assert iter.reader.index !== nothing
    chunks = GenomicFeatures.Indexes.overlapchunks(iter.reader.index.index, refindex, iter.interval)
    if !isempty(chunks)
        seek(iter.reader, first(chunks).start)
    end
    state = BAM.OverlapIteratorState(refindex, chunks, 1, BAM.Record())
    return iterate(iter, state)
end

function Base.iterate(iter::OverlapIterator, state)
    while state.chunkid ≤ lastindex(state.chunks)
        chunk = state.chunks[state.chunkid]
        while BGZFStreams.virtualoffset(iter.reader.stream) < chunk.stop
            read!(iter.reader, state.record)
            aln_interval = get_alignment_interval(state.record, iter.max_frag_size, iter.strand_type)
            aln_interval === nothing && continue
            c = GenomicFeatures.compare_overlap(aln_interval, Interval(iter.refname, iter.interval), isless)
            if c == 0
                return copy(state.record), state
            end
            # Default Base.iterate function from BioAlignments.BAM.overlaps.jl breaks
            # once the record is completely to the right of the GenomicFeatures interval,
            # but this causes some fragments to be missed, where the record is on the
            # negative strand and extends to the left into the GenomicFeatures interval
        end
        state.chunkid += 1
        if state.chunkid ≤ lastindex(state.chunks)
            seek(iter.reader, state.chunks[state.chunkid].start)
        end
    end
    return nothing
end
