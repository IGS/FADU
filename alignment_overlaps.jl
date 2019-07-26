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
    pp_only::Bool
    strand_type::String
    max_frag_size::UInt
end

function GenomicFeatures.eachoverlap(reader::BAM.Reader, interval::Interval, pp_only::Bool=false, strand_type::String="no", max_frag_size::UInt=1000)
    return GenomicFeatures.eachoverlap(reader, interval.seqname, interval.first:interval.last, pp_only, strand_type, max_frag_size)
end

function GenomicFeatures.eachoverlap(reader::BAM.Reader, interval, pp_only::Bool=false, strand_type::String="no", max_frag_size::UInt=1000)
    return GenomicFeatures.eachoverlap(reader, convert(Interval, interval), pp_only, strand_type, max_frag_size)
end

function GenomicFeatures.eachoverlap(reader::BAM.Reader, refname::AbstractString, interval::UnitRange, pp_only::Bool=false, strand_type::String="no", max_frag_size::UInt=1000)
    return OverlapIterator(reader, String(refname), interval, pp_only, strand_type, max_frag_size)
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
            record_type = determine_record_type(state.record, iter.max_frag_size, iter.pp_only)
            record_type === nothing && continue  # Skip alignments not being considered
            aln_interval = get_alignment_interval(state.record, record_type, is_reverse_stranded(iter.strand_type))
            c = GenomicFeatures.compare_overlap(aln_interval, Interval(iter.refname, iter.interval), isless)
            if c == 0
                #return aln_interval, state
                return copy(state.record), state
            end
        end
        state.chunkid += 1
        if state.chunkid ≤ lastindex(state.chunks)
            seek(iter.reader, state.chunks[state.chunkid].start)
        end
    end
    return nothing
end
