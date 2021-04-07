module Junctions

using ..Blocks
using ..Nodes
using ..Graphs: Counter, add!

export Junction
export junctions, left, right

struct Junction
    left::Node{Block}
    right::Node{Block}
end

blocks(j::Junction) = (left=j.left.block, right=j.right.block)

# --------------------------------
# extension of base operators

Base.hash(j::Junction) = Base.hash(Set([blocks(j), blocks(reverse(j))]))
Base.isequal(j₁::Junction, j₂::Junction) = Base.isequal(blocks(j₁), blocks(j₂)) || 
                                           Base.isequal(blocks(j₁), blocks(reverse(j₂)))

# --------------------------------
# custom operators

function left(j::Junction)
    return (
        block  = j.left.block,
        strand = j.left.strand,
    )
end

function right(j::Junction)
    return (
        block  = j.right.block,
        strand = j.right.strand,
    )
end

function reverse(j::Junction)
    return Junction(
        Node{Block}(j.right.block, !j.right.strand),
        Node{Block}(j.left.block, !j.left.strand),
    )
end

function junctions(paths)
    js = Dict{Junction, Counter}() 
    for path in paths
        for (i, node) ∈ enumerate(path.node[2:end])
            key = Junction(path.node[i], node)
            if key ∈ keys(js)
                add!(js[key], path.name)
            else
                js[key] = Counter([path.name => 1])
            end
        end
    end

    return js
end

end
