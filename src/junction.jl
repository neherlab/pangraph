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

blocks(j::Junction) = (left=(j.left.block,j.left.strand), right=(j.right.block,j.right.strand))

# --------------------------------
# extension of base operators

Base.hash(j::Junction) = Base.hash(blocks(j)) ⊻ Base.hash(blocks(reverse(j)))
Base.isequal(j₁::Junction, j₂::Junction) = Base.isequal(blocks(j₁), blocks(j₂)) || Base.isequal(blocks(j₁), blocks(reverse(j₂)))

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
        Node{Block}(j.left.block,  !j.left.strand),
    )
end

function junctions(paths)
    js = Dict{Junction, Counter}() 
    function push(key, name)
        if key ∈ keys(js)
            add!(js[key], name)
        else
            js[key] = Counter([name => 1])
        end
    end

    for path in paths
        for (i, node) ∈ enumerate(path.node[2:end])
            push(Junction(path.node[i], node), path.name)
        end

        if path.circular && length(path) > 1
            push(Junction(path.node[end],path.node[1]), path.name)
        end
    end

    return js
end

end
