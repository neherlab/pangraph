module Blocks

using Rematch, FStrings

using ..Utility: random_id, uncigar, wcpair, homologous, Alignment
using ..Nodes

import ..Graphs: pair, reverse_complement

export SNPMap, IndelMap #aux types
export Block 
export sequence, combine, swap! #operators

# ------------------------------------------------------------------------
# Block data structure

# aliases
SNPMap   = Dict{Int,UInt8}
IndelMap = Dict{Int,Union{Array{UInt8},Int}} # array = insertion; integer = deletion

mutable struct Block
    uuid::String
    sequence::Array{UInt8}
    mutation::Dict{Node{Block},SNPMap}
    indel::Dict{Node{Block},IndelMap}
end

# ---------------------------
# constructors

# simple helpers
Block(sequence,mutation,indel) = Block(random_id(),sequence,mutation,indel)
Block(sequence)                = Block(sequence, Dict{Node{Block},SNPMap}(), Dict{Node{Block},IndelMap}())
Block()                        = Block(UInt8[])

translate(dict, δ) = Dict(key=>Dict(x+δ => v for (x,v) in val) for (key,val) in dict)
function translate!(dict, δ)
    for (key, val) in dict
        dict[key] = Dict(x+δ => v for (x,v) in val)
    end
end

# TODO: rename to concatenate?
# serial concatenate list of blocks
function Block(bs::Block...)
    @assert all([isolates(bs[1]) == isolates(b) for b in bs[2:end]])

    sequence = join([b.sequence for b in bs])
    mutation = bs[1].mutation
    indel    = bs[1].indel

    δ = length(bs[1])
    for b in bs[2:end]
        merge!(mutation, translate(b.mutation, δ))
        merge!(indel, translate(b.indel, δ))

        δ += length(b)
    end

    return Block(sequence,mutation,indel)
end

# TODO: rename to slice?
# returns a subslice of block b
function Block(b::Block, slice)
    @show slice
    @assert slice.start >= 1 && slice.stop <= length(b)
    sequence = b.sequence[slice]
    mutation = translate(b.mutation, -slice.start)
    indel    = translate(b.indel, -slice.start)

    return Block(sequence,mutation,indel)
end

# ---------------------------
# operations

# simple operations
depth(b::Block) = length(b.mutation)
pair(b::Block)  = b.uuid => b

Base.length(b::Block) = Base.length(b.sequence)
function Base.length(b::Block, n::Node)
    length(x::Int)          = -x         # deletion
    length(x::Array{UInt8}) = length(x)  # insertion

    return length(b) + sum(length(v) for v in b.indel[n])
end

Base.show(io::IO, b::Block) = Base.show(io, (id=b.uuid, depth=depth(b)))

# complex operations
function reverse_complement(b::Block)
    seq = reverse_complement(b.sequence)
    len = length(seq)

    revcmpl(seq::Array{UInt8}) = reverse_complement(seq)
    revcmpl(n::Int)            = n
    revcmpl(dict::SNPMap)      = Dict(len-locus+1:wcpair[nuc] for (locus,nuc) in dict)
    revcmpl(dict::IndelMap)    = Dict(len-locus+1:revcmpl(val) for (locus,revcmpl) in dict)

    mutation = Dict(node => revcmpl(snps)  for (node, snps) in b.mutation)
    indel    = Dict(node => revcmpl(indel) for (node, indel) in b.indel)

    return Block(seq,mutation,indel)
end

# returns the sequence WITH mutations and indels applied to the consensus for a given tag 
function sequence(b::Block, node::Node{Block}; gaps=false)
    seq = copy(b.sequence)
    for (locus, snp) in b.mutation[node]
        seq[locus] = snp
    end

    indel(seq, locus, insert::Array{UInt8}) = cat(seq[1:locus], insert, seq[locus+1:end])
    indel(seq, locus, delete::Int)          = (gaps ? cat(seq[1:locus], '-'^delete, seq[locus+delete:end])
                                                    : cat(seq[1:locus], seq[locus+delete:end]))

    for locus in sort(keys(b.indel[node]),rev=true)
        indel = b.indel[node][locus]
        seq   = indel(seq, locus, indel)
    end

    return seq
end

function Base.append!(b::Block, node::Node{Block}, snp::SNPMap, indel::IndelMap)
    @assert node ∉ keys(b.mutation)
    @assert node ∉ keys(b.indel)

    b.mutation[node] = snp
    b.indel[node]    = indel
end

function swap!(b::Block, oldkey::Node{Block}, newkey::Node{Block})
    b.mutation[newkey] = pop!(b.mutation, oldkey)
    b.indel[newkey]    = pop!(b.indel, oldkey)
end

function swap!(b::Block, oldkey::Array{Node{Block}}, newkey::Node{Block})
    mutation = pop!(b.mutation, oldkey[1])
    indel    = pop!(b.indel, oldkey[1])

    for key in oldkey[2:end]
        merge!(mutation, pop!(b.mutation, key))
        merge!(indel, pop!(b.indel, key))
    end

    b.mutation[newkey] = mutation
    b.indel[newkey]    = indel
end

function reconsensus!(b::Block)
    aln = zeros(UInt8, depth(b), maximum(length(b, n) for n in keys(b.indel)))
end

function combine(qry::Block, ref::Block, aln::Alignment; maxgap=500)
    sequences,intervals,mutations,indels = homologous(
                                                uncigar(aln.cigar),
                                                qry.sequence,
                                                ref.sequence,
                                                maxgap=maxgap
                                           )

    blocks = NamedTuple{(:block,:kind),Tuple{Block,Symbol}}[]
    for (seq,pos,snp,indel) in zip(sequences,intervals,mutations,indels)
        @match (pos.qry, pos.ref) begin
            ( nothing, rₓ )  => push!(blocks, (block=Block(ref, rₓ), kind=:ref))
            ( qₓ , nothing ) => push!(blocks, (block=Block(qry, qₓ), kind=:qry))
            ( qₓ , rₓ )      => begin
                @assert !isnothing(snp)
                @assert !isnothing(indel)

                # slice both blocks
                r = Block(ref, rₓ)
                q = Block(qry, qₓ)

                # apply global snp and indels to all query sequences
                # XXX: do we have to worry about overlapping insertions/deletions?
                for node in keys(q.mutation)
                    merge!(q.mutation[node], snp)
                    merge!(q.indel[node], indel)
                end

                new = Block(seq,snp,indel)
                reconsensus!(new)

                push!(blocks, (block=new, kind=:all))
            end
        end
    end

    return blocks
end

end
