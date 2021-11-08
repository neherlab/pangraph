module Edges

using ..Paths
using ..Nodes
using ..Blocks
using ..Graphs

export edges

struct Position
    path  :: Path
    node  :: Tuple{Node{Block},Node{Block}}
    index :: Tuple{Int,Int} # positions on path
    locus :: Int # breakpoint on sequence
end

mutable struct Edge
    block  :: Tuple{Block, Block}
    invert :: Bool # changes strand
    nodes  :: Array{Position}
end

struct EdgeSet
    edges :: Array{Edge,1}
    links :: Dict{Node{Block}, Tuple{Maybe{Edge}, Maybe{Edge}}}
end

function equals(e::Edge, block1::Block, block2::Block, invert::Bool)
    return e.block[1] == block1 && e.block[2] == block2 && invert == e.invert
end

Edge(block1::Block, block2::Block, invert::Bool, pos::Position) = Edge((block1, block2), invert, Position[pos])

function edges(G)
    edgeset = Edge[]

    # XXX: think of a more declarative syntax than a raw tuple
    function addedge!(node, link, path, locus)
        if node[1].block.uuid > link[1].block.uuid
            node, link = link, node
        end

        invert   = node[1].strand ⊻ link[1].strand
        position = Position(
            path,
            (node[1],link[1]),
            (node[2],link[2]),
            locus
        )
        for edge in edgeset
            if equals(edge, node[1].block, link[1].block, invert)
                push!(edge.nodes, position)
                return
            end
        end
        # if we are here, then this is a new edge
        push!(edgeset, Edge(node[1].block, link[1].block, invert, position))
    end

    for path in values(G.sequence)
        for (i,link) in enumerate(path.node[2:end])
            node = path.node[i]
            addedge!((node,i), (link,i+1), path, path.position[i])
        end

        if path.circular
            addedge!((path.node[end], length(path.node)), (path.node[1],1), path, path.position[end])
        end
    end

    return edgeset
end

end
