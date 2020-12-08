module Graphs

using Match, FStrings
using GZip # NOTE: for debugging purposes

import JSON

# ------------------------------------------------------------------------
# imports w/ exported global dispatch preamble

export pair
function pair(item) end

include("counter.jl")
include("util.jl")
# NOTE: commented out during debugging stage
#       will move to a module file later
# include("pool.jl")
include("node.jl")
include("block.jl")
include("path.jl")
include("junction.jl")

using .Utility: read_fasta, name, columns
using .Nodes
using .Blocks
using .Paths
using .Junctions

export Graph
export graphs, marshal, serialize

# ------------------------------------------------------------------------
# graph data structure

struct Graph
    block::Dict{String,Block}   # uuid      -> block
    sequence::Dict{String,Path} # isolation -> path
    # TODO: add edge/junction data structure?
end

include("align.jl")
using .Align

# --------------------------------
# constructors

function Graph(name::String, sequence::Array{UInt8}; circular=false)
    block = Block(sequence)
    path  = Path(name, Node{Block}(block); circular=circular)

    add!(block, path.node[1], SNPMap(), IndelMap())

    return Graph(
         Dict([pair(block)]), 
         Dict([pair(path)]),
         # TODO: more items...
    )
end

graphs(io::IO) = [Graph(name(record), record.seq) for record in read_fasta(io)]

# --------------------------------
# operators

# XXX: break into smaller functions
#      too long
function detransitive!(G::Graph)
    isosᵥ = count_isolates(values(G.sequence))

    # collect all junctions that transitively pass isolates through
    transitives = Junction[]
    for (j, isosⱼ) in junctions(values(G.sequence))
        if ((isosᵥ[j.left.block] == isosᵥ[j.right.block]) &&
            (isosᵥ[j.left.blockj] == isosⱼ))
            push!(transitives, j)
        end
    end

    Link  = NamedTuple{(:block,:strand),(Block, Bool)}
    Chain = Array{Link}
    
    rev(l::Link)  = (block=l.block,strand=!l.strand)
    rev(c::Chain) = [rev(b) for b in reverse(c)]

    # build chains by threading consecutive transitive junctions
    # TODO: audit each line carefully
    chain = Dict{Block, Chain}()
    for j in transitives
        if j.left.block ∈ keys(chain) && j.right.block ∈ keys(chain)
            c₁, c₂ = chain[j.left.block], chain[j.right.block]
            if c1 == c2
                continue
            end

            newc = Block[]
            if left(j) == c₁[end] && right(j) == c₂[1]
                newc = cat(c₁, c₂, dims=1)
            elseif left(j) == c₁[end] && rev(right(j)) == c₂[end]
                newc = cat(c₁, rev(c₂), dims=1)
            elseif rev(left(j)) == c₁[1] && right(j) == c₂[1]
                newc = cat(rev(c₁), c₂, dims=1)
            elseif rev(left(j)) == c₁[1] && rev(right(j)) == c₂[end]
                newc = cat(c₂, c₁, dims=1)
            else
                error("case not covered")
            end

            for b in newc
                chain[b] = newc
            end

        elseif j.left.block ∈ keys(chain)
            c₀ = chain[j.left.block]
            if left(j) == c₀[end]
                push!(c₀, right(j))
            elseif rev(j.left.block) == c₀[1]
                pushfirst!(c₀, rev(right(j)))
            else
                error("chains should be linear")
            end
            chain[j.right.block] = c₀

        elseif j.right.block ∈ keys(chain)
            c₀ = chain[j.right.block]
            if right(j) == c₀[end]
                push!(c₀, rev(left(j)))
            elseif right(j) == c₀[1]
                pushfirst!(c₀, left(j))
            else
                error("chains should be linear")
            end
            chain[j.left.block] = c₀

        else
            chain[j.left.block]  = [left(j), right(j)]
            chain[j.right.block] = chain[j.left.block] 
        end
    end

    # merge chains into one block
    for c in Set(values(chain))
        isos = isosᵥ[c[1].block]
        new  = Block([s ? b : reverse_complement(b) for (b,s) in c]...)

        for iso in isos
            replace!(G.sequence[iso], c, new)
        end

        for (b,_) in c
            pop!(G.blocks, b)
        end
    end
end

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

function marshal_json(io, G::Graph)
    out = {}
end

function marshal(io, G::Graph; fmt=:fasta)
    @match fmt begin
        :fasta || :fa => marshal_fasta(io, G)
        :json         => marshal_json(io, G)
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
        println(">aligning...")
        graph    = align(isolates[1], isolates[2])
    end
end

end
