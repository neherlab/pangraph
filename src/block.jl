module Blocks

using Match, FStrings

using ..Utility: random_id
using ..Nodes

import ..Graphs: pair

export SNPMap, IndelMap #aux types
export Block 
export sequence, merge, add! #operators

# ------------------------------------------------------------------------
# Block data structure

# aliases
SNPMap   = Dict{Int,UInt8}
IndelMap = Dict{Int,Union{Array{UInt8},Int}} # array = insertion; integer = deletion

struct Block
    uuid::String
    sequence::Array{UInt8}
    mutation::Dict{Node{Block},SNPMap}
    indel::Dict{Node{Block},IndelMap}
end

# ---------------------------
# constructors

# simple helpers
Block() = Block(random_id(), UInt8[], Dict{Node{Block},SNPMap}(), Dict{Node{Block},IndelMap}())
Block(sequence) = Block(random_id(), sequence, Dict{Node{Block},SNPMap}(), Dict{Node{Block},IndelMap}())
Block(sequence,mutation,indel) = Block(random_id(),sequence,mutation,indel)

# serial concatenate list of blocks
function Block(bs::Block...)
    @assert all([isolates(bs[1]) == isolates(b) for b in bs[2:end]])

    sequence = join([b.sequence for b in bs])
    mutation = bs[1].mutation
    indel    = bs[1].insertion

    # TODO: finish the concatenation of mutations/etc dictionaries
    #       need to thread degeneracies through the list correctly
    δ = length(bs[1])
    for b in bs[2:end]
        δ += length(b)
    end

    return Block(sequence,mutation,indel)
end

# returns a subslice of block b
function Block(b::Block, range::UnitRange{Int})
    @assert range.start >= 1 && range.stop <= length(b)
    translate(dict) = Dict(iso => Dict(x-range.start => val for (x,val) in d) for (iso,d) in dict)
    sequence = b.sequence[range]
    mutation = translate(b.mutation)
    indel    = translate(b.indel)

    return Block(sequence,mutation,indel)
end

# ---------------------------
# operations

# simple operations
length(b::Block) = length(b.sequence)
depth(b::Block)  = length(b.mutation)
pair(b::Block)   = b.uuid => b

# complex operations

# NOTE: this is deprecated now that we switched from tags to nodes as keys!
#       when do we use this? do we need it anymore since the switch?
#
# # returns the count of isolates within the block
# function isolates(b::Block)
#     count = Dict{String,Int}()
#     for node in b.mutation
#         if node.isolate in keys(count)
#             count[tag.isolate] += 1
#         else
#             count[tag.isolate]  = 1
#         end
#     end

#     return count
# end

# returns the sequence WITH mutations and indels applied to the consensus for a given tag 
function sequence(b::Block, node::Node{Block}; gaps=false)
    seq = copy(b.sequence)
    # mutations
    for (locus, snp) in b.mutation[node]
        seq[locus] = snp
    end

    # indels
    indel(seq, locus, insert::Array{UInt8}) = cat(seq[1:locus], insert, seq[locus+1:end])
    indel(seq, locus, delete::Int) = (gaps ? cat(seq[1:locus], '-'^delete, seq[locus+delete:end])
                                           : cat(seq[1:locus], seq[locus+delete:end]))

    for locus in sort(keys(b.indel[node]),rev=true)
        indel = b.indel[node][locus]
        seq   = indel(seq, locus, indel)
    end

    return seq
end

function add!(b::Block, node::Node{Block}, snp::SNPMap, indel::IndelMap)
    @assert node ∉ keys(b.mutation)
    @assert node ∉ keys(b.indel)
    b.mutation[node] = snp
    b.indel[node]    = indel
end

end
