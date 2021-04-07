module Paths

using ..Nodes
using ..Blocks

import ..Graphs: 
    pair,
    Counter, add!

export Path
export count_isolates

struct Path
    name     :: String
    node     :: Array{Node{Block}}
    offset   :: Union{Int,Nothing}
    circular :: Bool
end

# --------------------------------
# constructors

Path(name::String,node::Node{Block};circular::Bool=false) = Path(name,[node],circular ? 0 : nothing,circular)

# --------------------------------
# operators

pair(p::Path) = p.name => p
Base.show(io::IO, p::Path) = Base.show(io, (name=p.name, blocks=p.node))

# used for block merging
function Base.replace!(p::Path, old::Block, new::Array{Block})
    indices = Int[]
    inserts = Array{Node{Block}}[]
    for (i, n₁) in enumerate(p.node)
        if n₁.block != old
            continue
        end

        push!(indices, i)

        nodes = (n₁.strand ? [Node{Block}(nb) for nb in new] 
                           : [Node{Block}(nb;strand=false) for nb in reverse(new)])
        for (blk,n₂) in zip(n₁.strand ? new : reverse(new), nodes)
            swap!(blk, n₁, n₂)
        end

        push!(inserts, nodes)
    end

    # reverse so that lower indices don't shift upper while we iterate
    reverse!(indices)
    reverse!(inserts)

    for (index, nodes) in zip(indices, inserts)
        splice!(p.node, index, nodes)
    end
end

# XXX: should we create an interval data structure that unifies both cases?
# XXX: wrap as an iterator instead of storing the whole array in memory?
# XXX: do we have to deal with reverse chain case?
intervals(starts, stops) = [start:stop for (start,stop) in zip(starts,stops)]

function intervals(starts, stops, gap, len)
    if (stops[1]-starts[1]) == gap
        return intervals(starts, stops)
    elseif (len-starts[end]+stops[1]+1) == gap
        return [(starts[end]:len; 1:stops[1]);[start:stop for (start,stop) in zip(starts[1:end-1], stops[2:end])]]
    else
        error("unrecognized gap pattern in block replacement of path") 
    end
end

# used for detransitive
const Link = NamedTuple{(:block, :strand), Tuple{Block, Bool}}
function Base.replace!(p::Path, old::Array{Link}, new::Block)
    start, stop = old[1].block, old[end].block

    i₁ₛ = findall((n) -> n.block == start, p.node)
    i₂ₛ = findall((n) -> n.block == stop,  p.node)

    @assert length(i₁ₛ) == length(i₂ₛ)

    # XXX: the ordering here is wrong for circular genomes...
    iₛ = p.circular ? intervals(i₁ₛ, i₂ₛ, length(old), length(p.node)) : intervals(i₁ₛ, i₂ₛ)
    sₛ = map(i₁ₛ) do i
        p.node[i].strand == old[1].strand
    end

    splice!(nodes, i::AbstractArray, new)                       = Base.splice!(nodes, i, new)
    splice!(nodes, i::Tuple{AbstractArray,AbstractArray}, new)  = let 
        Base.splice!(nodes, i[1], new)
        Base.splice!(nodes, i[2])
    end

    # NOTE: have to correct for the overhang in the case of circular intervals
    swap!(b, i::AbstractArray, new)                       = Blocks.swap!(b, p.node[i], new)
    swap!(b, i::Tuple{AbstractArray,AbstractArray}, new)  = Blocks.swap!(b, p.node[[length(p.node)-length(i[1]):length(p.node);i[2]]], new)

    for (i,s) ∈ zip(reverse(iₛ), reverse(sₛ))
        node = Node(new, s)

        swap!(new, i, node)
        splice!(p.node, i, [node])
    end
end

# XXX: store as field in block?
#      conversely we could store a pointer to the Path in each node.
#      this would allow us to access which isolate each node belongs to.
#      for now we keep it global to have memory usage be minimal.
#      important to consider
# XXX: strongly type iterator for more specificity?
function count_isolates(paths)
    blocks  = Dict{Block, Counter}()

    for path in paths
        for node in path.node
            if node.block ∈ keys(blocks)
                add!(blocks[node.block], path.name)
            else
                blocks[node.block] = Counter([path.name => 1])
            end
        end
    end

    return blocks
end

end
