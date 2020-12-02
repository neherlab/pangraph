module Pangraph

using Match, FStrings

# NOTE: commented out during debugging stage
include("util.jl")
# include("pool.jl")
# include("align.jl")
include("block.jl")
include("path.jl")

using .Utility: read_fasta, name
using .Blocks, .Paths

export Graph, Graphs, marshal, serialize

# ------------------------------------------------------------------------
# graph data structure

struct Graph
    block::Dict{String,Block}
    sequence::Dict{String,Path}
    # TODO: add edge data structure
end

# --------------------------------
# constructors

function Graph(name::String, sequence::Array{Char})
    block = Block(name, sequence)
    path  = Path(name, Node(block))

    return Graph(
         Dict([(block.uuid, block)]), 
         Dict([(path.name, path)])
    )
end

Graphs(io::IO) = [Graph(name(record), record.seq) for record in read_fasta(io)]

# --------------------------------
# operators

# TODO: think of better names
Nₛ(G::Graph) = length(G.sequence)
Nᵥ(G::Graph) = length(G.block)

# ------------------------------------------------------------------------
# serialization

# helper functions w/ common functionality
function columns(s; nc=80)
    nr   = ceil(Int64, length(s)/nc)
    l(i) = 1+(nc*(i-1)) 
    r(i) = min(nc*i, length(s))
    rows = [String(s[l(i):r(i)]) for i in 1:nr]
    return join(rows,'\n')
end

function write_fasta(io, name, seq)
    write(io, '>', name, '\n')
    write(io, columns(seq), '\n')
end

function marshal_fasta(io, G::Graph)
    for (i, b) in enumerate(values(G.block))
        write_fasta(io, b.uuid, b.sequence)
    end
end

function marshal(io, G::Graph; fmt=:fasta)
    @match fmt begin
        :fasta || :fa => marshal_fasta(io, G)
        _ => error(f"{format} not a recognized output format")
    end
end

# TODO: can we generalize to multiple individuals
#       equivalent to "highway" detection
function serialize(io, G::Graph)
    if length(G.sequence) != 1
        error("only singleton graphs implemented")
    end

    name = collect(keys(G.sequence))[1]
    seq  = collect(values(G.block))[1].sequence

    write_fasta(io, name, seq)
end

end
