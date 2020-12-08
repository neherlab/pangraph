module Paths

using ..Nodes
using ..Blocks

import ..Graphs: pair, Counter, add!

export Path
export count_isolates

struct Path
    name::String
    node::Array{Node{Block}}
    offset::Union{Int,Nothing}
    circular::Bool
end

# --------------------------------
# constructors

Path(name::String,node::Node{Block};circular::Bool=false) = Path(name,[node],circular ? 0 : nothing,circular)

# --------------------------------
# operators

pair(p::Path) = p.name => p

# XXX: do we just want to build a brand new array incrementally?
function Base.replace!(p::Path, old::Block, new::Array{Block})
    indices = Int[]
    inserts = Array{Node{Block}}[]
    for (i, n) in enumerate(p.node)
        if n.block != old
            continue
        end

        push!(indices, i)

        nodes = (n.strand ? [Node{Block}(nb) for nb in new] 
                          : [Node{Block}(nb;strand=false) for nb in reverse(new)])
        for new in nodes
            swap!(old, n, new)
        end

        push!(inserts, nodes)
    end

    reverse!(indices)
    reverse!(inserts)

    for (index, nodes) in zip(indices, inserts)
        splice!(p.node, index, nodes)
    end
end

# XXX: should we create an interval data structure that unifies both cases?
# XXX: wrap as an iterator instead of storing the whole array in memory?
intervals(starts, stops) = [stop:-1:start for (start,stop) in zip(starts,stops)]

function intervals(starts, stops, gap, len)
    if (stops[1]-starts[1]) == gap
        return intervals(starts, stops)
    elseif (len-starts[end]+stops[1]+1) == gap
        return [(starts[end]:len; 1:stops[1]);[start:stop for (start,stop) in zip(starts[1:end-1], stops[2:end])]]
    else
        error("unrecognized gap pattern in block replacement of path") 
    end
end

Link = NamedTuple{(:block, :strand), Tuple{Block, Bool}}
function Base.replace!(p::Path, old::Array{Link}, new::Block)
    start, stop = old[1].block, old[end].block

    i₁ₛ = findall((n) -> n.block == start, p.node)
    i₂ₛ = findall((n) -> n.block == stop,  p.node)

    @assert length(i₁ₛ) == length(i₂ₛ)

    s  = new.strand == old[1].strand
    iₛ = (p.circular ? intervals(i₁ₛ, i₂ₛ, length(old), length(p.node))
                     : intervals(i₁ₛ, i₂ₛ))

    Interval         = UnitRange{Int}

    splice!(nodes::Array{Node{Block}}, i::UnitRange{Int}, new::Node{Block}) = Base.splice!(nodes, i, new)
    splice!(nodes::Array{Node{Block}}, i::Tuple{UnitRange{Int},UnitRange{Int}}, new::Node{Block}) = Base.splice!(nodes, i, [1], new); Base.splice!(nodes, idx[2])

    # NOTE: have to correct for the overhang in the case of circular intervals
    swap!(b::Block, i::UnitRange{Int}, new::Node{Block}) = Block.swap!(b, p.nodes[i], new)
    swap!(b::Block, i::Tuple{UnitRange{Int},UnitRange{Int}}, new::Node{Block}) = Block.swap!(b, p.nodes[[length(p.nodes)-length(i[1]):length(p.nodes);i[2]]], new)

    for interval in reverse(iₛ)
        n = Node(new, s)
        splice!(p.nodes, interval, n)
        swap!(new, interval, n)
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
