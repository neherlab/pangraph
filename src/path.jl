module Paths

using ..Nodes
using ..Blocks

import ..Graphs: pair, Counter, add!

export Path
export replace!

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
function replace!(p::Path, old::Block, new::Array{Block})
    indices = Int[]
    inserts = Array{Node}[]
    for (i, n) in p.node
        if n.block != old
            continue
        end

        push!(indices, i)

        nodes = (n.strand ? [Node(nb) for nb in new] 
                          : [Node(nb;strand=false) for nb in reverse(new)])
        for new in nodes
            swap!(m.block, n, new)
        end

        push!(inserts, nodes)
    end

    for (index,nodes) in reverse(zip(indices, inserts))
        splice!(p.nodes, index, nodes)
    end
end

Link = NamedTuple{(:block, :strand), (Block, Bool)}
function replace!(p::Path, old::Array{Link}, new::Block)
    start, stop, len = old[1].block, old[end].block, length(old)

    i₁ₛ = findall((n) -> n.block == start, p.node)
    i₂ₛ = findall((n) -> n.block == stop,  p.node)

    # XXX: pair the i₁ & i₂ into (i₁, i₂) s.t. |i₂-i₁| = len
    #      take into consideration the circularity of path p
    #      deal with strand sign correctly here

    # XXX: for (i₁, i₂) in pairs
    #           splice!(p.nodes, i₁:i₂, Node(new)) :: Take into consideration circularity!
    #           swap out p.nodes[i₁:i₂] in Block mutation & indel dictionaries for newly created node
    #      end

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
