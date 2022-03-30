module GFA

import Base: print
import ..Graphs: Graph, sequence, marshal_gfa

"""
    struct Segment
        name     :: String
        sequence :: Array{UInt8}
        depth    :: Int
    end

Store a GFA segment, i.e. an edge of an alignment graph that holds a contiguous sequence
Depth is the number of genomes, including duplications, that contain the sequence.
"""
struct Segment
    name     :: String
    sequence :: Array{UInt8}
    depth    :: Int
end

function print(io::IO, s::Segment)
    print(io, "S\t$(s.name)\t*\tLN:i:$(length(s.sequence))\tRC:i:$((s.depth*length(s.sequence)))")
end

const Node = NamedTuple{
    (:segment,:orientation),
    Tuple{Segment,Bool}
}

struct Key
    prev :: Node
    curr :: Node
end

Base.hash(k::Key) = Base.hash(Set(
    [(prev=k.prev,curr=k.curr),(prev=(segment=k.curr.segment,orientation=!k.curr.orientation),curr=(segment=k.prev.segment,orientation=!k.prev.orientation))]
))
# Base.hash(k::Key, h::UInt) = Base.hash(Base.hash(k), h)

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

"""
    struct Path
        name     :: String
        segments :: Array{Node,1}
        circular :: Bool
    end

Store a GFA path, i.e. a sequence of segments that represents an observed genome.
"""
struct Path
    name     :: String
    segments :: Array{Node,1}
    circular :: Bool
end

function print(io::IO, nodes::Array{Node,1})
    print(io, join(("$(n.segment.name)$(polarity(n.orientation))" for n in nodes ), ","))
end

circular(flag::Bool) =  flag ? "TP:Z:circular" : "TP:Z:linear"

function print(io::IO, p::Path)
    print(io, "P\t$(p.name)\t$(p.segments)\t$(circular(p.circular))*")
end

"""
    marshal_gfa(io::IO, G::Graph; opt=nothing)

Output pangraph `G` to IO stream `io`.
`opt` can include two functions, to be accessed in fields `connect` and `output`.
`connect` is a function that takes a node and returns true or false if it should be connected in the GFA output.
`output` is an equivalent function signature, but controls whether the node is output at all.
"""
function marshal_gfa(io::IO, G::Graph; opt=nothing)
    # unpack options
    filter = (
        opt === nothing
        ? (
            connect = (node) -> false,
            output  = (segment, isolates) -> false,
          )
        : (
            connect = opt.connect,
            output  = opt.output,
          )
    )

    # wrangle data
    segments = Dict(
        let
            segment = Segment(
                block.uuid,
                sequence(block),
                length(block.mutate)
            )
            isolates = [
                let
                    isolate = ""
                    for (name, path) in G.sequence
                        if key ∈ path.node
                            isolate = name
                            break
                        end
                    end
                    isolate
                end for key in keys(block.mutate)
            ]
            block => filter.output(segment,isolates) ? nothing : segment
        end for block in values(G.block)
    )

    links = Dict{Key, Link}()
    paths = Array{Path,1}(undef, length(G.sequence))

    addlink! = (prev, curr) -> let
        key = Key(prev, curr)
        if key ∉ keys(links)
            links[key] = Link(prev, curr, 1)
        else
            links[key].depth += 1
        end
    end

    for (i,path) in enumerate(values(G.sequence))
        nodes = [
            let
                segment = segments[node.block]
                if segment !== nothing && !filter.connect(node)
                    (
                     segment     = segment,
                     orientation = node.strand,
                    )
                else
                    nothing
                end
            end
            for node in path.node
        ]
        filter!(n->n!==nothing, nodes)

        if length(nodes) == 0
            continue
        end

        for (j,node) in enumerate(nodes[2:end])
            addlink!(nodes[j], node)
        end

        if path.circular
            addlink!(nodes[end], nodes[1])
        end

        paths[i] = Path(path.name, nodes, path.circular)
    end

    # export data
    write(io, "H\tVN:Z:1.0\n")

    write(io, "# pancontigs\n")
    for segment in values(segments)
        segment === nothing && continue

        write(io, "$(segment)\n")
    end

    write(io, "# edges\n")
    for link in values(links)
        write(io,"$(link)\n")
    end

    write(io, "# sequences\n")
    for i = 1:length(paths)
        if isassigned(paths, i)
            write(io, "$(paths[i])\n")
        end
    end
end

end
