module Edges

using ..Paths
using ..Nodes
using ..Blocks
using ..Graphs

using Infiltrator

export edges

struct Position
    path  :: Path
    node  :: Tuple{Node{Block},Node{Block}}
    index :: Tuple{Int,Int} # positions on path
    locus :: Int # breakpoint on sequence
end

Base.show(io::IO, x::Position) = Base.show(io,
    (
        path=x.path.name,
        blocks=(x.node[1].block, x.node[2].block),
        strand=(x.node[1].strand, x.node[2].strand),
        index=(x.index[1], x.index[2]),
        locus=x.locus
    )
)

mutable struct Edge
    block  :: Tuple{Block, Block}
    invert :: Bool # changes strand
    nodes  :: Array{Position}
end

Breakpoints = Array{Tuple{Int,Int},1}

function equals(e::Edge, block1::Block, block2::Block, invert::Bool)
    return e.block[1] == block1 && e.block[2] == block2 && invert == e.invert
end

Edge(block₁::Block, block₂::Block, invert::Bool, pos::Position) = Edge((block₁, block₂), invert, Position[pos])

#=
function connections(G::Graph, edges::Array{Edge,1})
    # collect all edges that flow into each vertex
end
=#

function edges(G)
    edgeset = Edge[]

    function addedge!(side₁, side₂, path, locus)
        if side₁.node.block.uuid > side₂.node.block.uuid
            side₁, side₂ = side₂, side₁
        end

        invert   = side₁.node.strand ⊻ side₂.node.strand
        position = Position(
            path,
            (side₁.node, side₂.node),
            (side₁.index,side₂.index),
            locus
        )
        for edge in edgeset
            if equals(edge, side₁.node.block, side₂.node.block, invert)
                push!(edge.nodes, position)
                return
            end
        end
        # if we are here, then this is a new edge
        push!(edgeset, Edge(side₁.node.block, side₂.node.block, invert, position))
    end

    for path in values(G.sequence)
        for (i,link) in enumerate(path.node[2:end])
            node = path.node[i]
            addedge!(
                (node=node,index=i),
                (node=link,index=i+1),
                path, path.position[i]
            )
        end

        if path.circular
            addedge!(
                (node=path.node[end],index=length(path.node)),
                (node=path.node[1],index=1),
                path, path.position[end]
            )
        end
    end

    return edgeset
end

function breakpoints(G)
    edgeset = edges(G)

    # collect all edges for each block, segregated by position
    vertex = Dict(
        block => (
            left  = Dict{Position,Edge}(),
            right = Dict{Position,Edge}()
        ) for block in values(G.block)
    )

    function connect(edge, block)
        for locus in edge.nodes
            i = locus.node[1] == block ? 1 : 2
            j = i == 1 ? 2 : 1
            N = length(locus.path.node)

            strand  = locus.node[i].strand
            toright = if (locus.index[j] % N) + 1 == locus.index[i]
                true
            elseif (locus.index[i] % N) + 1 == locus.index[j]
                false
            else
                # XXX: debugging assumptions
                # remove when confident
                error("bad index arithmetic")
            end

            if strand ⊻ toright
                vertex[block].left[locus] = edge
            else
                vertex[block].right[locus] = edge
            end
        end
    end

    for edge in edgeset
        block = edge.block
        connect(edge, block[1])
        connect(edge, block[2])
    end

    breakpoints = Dict{Tuple{String,String}, Breakpoints}()

    for block ∈ values(G.block) # iterate over all blocks
        for edge in vertex[block] # iterate over left & right breaks
            loci = collect(keys(edge))
            for (i,pos₁) in enumerate(loci) # all breakpoints on left/right edge
                e₁ = edge[pos₁]
                for pos₂ in loci[i+1:end]
                    e₂ = edge[pos₂]
                    if e₁ != e₂
                        key, val = if pos₁.path.name > pos₂.path.name
                            (pos₂.path.name, pos₁.path.name),
                            (pos₂.locus, pos₁.locus)
                        else
                            (pos₁.path.name, pos₂.path.name),
                            (pos₁.locus, pos₂.locus)
                        end
                        @infiltrate key == ("NZ_CP006918", "NZ_CP015822") && val[1] == 5242941
                        if key ∈ keys(breakpoints)
                            push!(breakpoints[key],val)
                        else
                            breakpoints[key] = [val]
                        end
                    end
                end
            end
        end
    end

    return breakpoints, function(name₁::AbstractString, name₂::AbstractString)
        name₁ == name₂ && return Breakpoints([])
        name₁ > name₂  && return breakpoints[(name₂,name₁)]
        return breakpoints[(name₁,name₂)]
    end
end

end
