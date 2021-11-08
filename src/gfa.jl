module GFA

import Base: print
import ..Graphs: Graph, sequence, marshal_gfa

struct Segment
    name     :: String
    sequence :: Array{UInt8}
    depth    :: Int
end

function print(io::IO, s::Segment)
    print(io, "S\t$(s.name)\t*\tLN:i:$(length(s.sequence))\tRC:i:$((s.depth*length(s.sequence)))")
end

Node = NamedTuple{
    (:segment, :orientation),
    Tuple{Segment,Bool}
}

polarity(x::Bool) = x ? "+" : "-"
function print(io::IO, n::Node)
    print(io, "$(n.segment.name)\t$(polarity(n.orientation))")
end

mutable struct Link
    from  :: Node
    to    :: Node
    depth :: Int
end

function print(io::IO, l::Link)
    print(io, "L\t$(l.from)\t$(l.to)\t*\tRC:i:$(l.depth)")
end

struct Path
    name     :: String
    segments :: Array{Node,1}
end

function print(io::IO, nodes::Array{Node,1})
    print(io, join(("$(n.segment.name)$(polarity(n.orientation))" for n in nodes ), ","))
end

function print(io::IO, p::Path)
    print(io, "P\t$(p.name)\t$(p.segments)\t*")
end

function marshal_gfa(io::IO, G::Graph)
    # wrangle data
    segments = Dict(
        block => Segment(
            block.uuid,
            sequence(block),
            length(block.mutate)
        ) for block in values(G.block)
    )

    links = Dict{Set{Node}, Link}()
    paths = Array{Path,1}(undef, length(G.sequence))

    addlink! = (prev, curr) -> let
        key = Set([prev, curr])
        if key ∉ keys(links)
            links[key] = Link(prev, curr, 1)
        else
            links[key].depth += 1
        end
    end

    for (i,path) in enumerate(values(G.sequence))
        nodes = [let
             if length(node.block.mutate) ≤ 30 & length(node.block.mutate) > 1
                (
                 segment     = segments[node.block],
                 orientation = node.strand,
                )
            else
                nothing
            end
        end for node in path.node]
        filter!(n->n!==nothing, nodes)

        for (j,node) in enumerate(nodes[2:end])
            addlink!(nodes[j], node)
        end

        if path.circular
            addlink!(nodes[end], nodes[1])
        end

        paths[i] = Path(path.name, nodes)
    end

    # export data
    write(io, "H\tVN:Z:1.0\n")

    write(io, "# pancontigs\n")
    for segment in values(segments)
        write(io, "$(segment)\n")
    end

    write(io, "# edges \n")
    for link in values(links)
        write(io,"$(link)\n")
    end

    write(io, "# sequences \n")
    for path in values(paths)
        write(io,"$(path)\n")
    end
end

end
