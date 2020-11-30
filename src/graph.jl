module Pangraph

using Match, FStrings

# NOTE: commented out during debugging stage
# include("util.jl")
# include("pool.jl")
# include("align.jl")
include("block.jl")
include("path.jl")

using .Blocks, .Paths

export Graph, write

# ------------------------------------------------------------------------
# graph data structure

struct Graph
    blocks::Dict{String,Block}
    sequence::Dict{String,Path}
    # TODO: add edge data structure
end

# --------------------------------
# constructors

function Graph(name::String, sequence::Array{Char})
    block = Block(name, sequence)
    path  = Path(name, Node(block))

    return Graph(
         Dict([(block.id, block)]), 
         Dict([(path.name, path)])
    )
end

# ------------------------------------------------------------------------
# serialization

function write_fasta(io, G::Graph; numcols=80)
    NL = '\n'
    for block in values(G.blocks)
        write(io, "f>{block.uuid}{NL}")
        write(io, join([block.sequence[1+numcols*i:numcols*(i+1)] for i in 1:ceil(length(block),numcols)], NL))
    end
end

function write(io, G::Graph; fmt=:fasta)
    @match fmt begin
        :fasta || :fa => write_fasta(io, G)
        _ => error(f"{format} not a recognized output format")
    end
end

end
