module Graphs

using Match, FStrings
using GZip # NOTE: for debugging purposes

export pair
function pair(item) end

# NOTE: commented out during debugging stage
#       will move to a module file later
include("util.jl")
# include("pool.jl")
# include("align.jl")
include("node.jl")
include("block.jl")
include("path.jl")

using .Utility: read_fasta, name, columns
using .Nodes
using .Blocks
using .Paths

export Graph
export graphs, marshal, serialize

# ------------------------------------------------------------------------
# graph data structure

struct Graph
    block::Dict{String,Block}   # uuid      -> block
    sequence::Dict{String,Path} # isolation -> path
    # TODO: add edge/junction data structure
end

# --------------------------------
# constructors

function Graph(name::String, sequence::Array{UInt8}; circular=false)
    println(">", name)
    block = Block(sequence)
    path  = Path(name, Node{Block}(block); circular=circular)

    println("> adding...")
    add!(block, path.node[1], SNPMap(), IndelMap())

    println("> building...")
    return Graph(
         Dict([pair(block)]), 
         Dict([pair(path)]),
         # TODO: more items...
    )
end

graphs(io::IO) = [Graph(name(record), record.seq) for record in read_fasta(io)]

# --------------------------------
# operators

# TODO: think of better names
Nₛ(G::Graph) = length(G.sequence)
Nᵥ(G::Graph) = length(G.block)

# ------------------------------------------------------------------------
# i/o & (de)serialization

# helper functions w/ common functionality
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

function test()
    GZip.open("data/generated/assemblies/isolates.fna.gz", "r") do io
        isolates = graphs(io)
        # graph    = align(isolates...)
    end
end

end
