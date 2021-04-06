module Graphs

using GZip # NOTE: for debugging purposes
using Random
using Rematch

import JSON

# ---------------------------
# functions to extend in submodules

export pair
function pair(item) end

export reverse_complement
function reverse_complement(item) end

include("interval.jl")
include("counter.jl")
include("util.jl")
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
export graphs, marshal, serialize, detransitive!

# ------------------------------------------------------------------------
# graph data structure

struct Graph
    block    :: Dict{String,Block}   # uuid      -> block
    sequence :: Dict{String,Path} # isolation -> path
    # TODO: add edge/junction data structure?
end

include("align.jl")
using .Align

# --------------------------------
# constructors

function Graph(name::String, sequence::Array{UInt8}; circular=false)
    block = Block(sequence)
    path  = Path(name, Node{Block}(block); circular=circular)

    append!(block, path.node[1], SNPMap(), InsMap(), DelMap())

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
Link  = NamedTuple{(:block,:strand),Tuple{Block, Bool}}
Chain = Array{Link}
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

    rev(l::Link)  = (block=l.block,strand=!l.strand)
    rev(c::Chain) = [rev(b) for b in reverse(c)]

    # build chains by threading consecutive transitive junctions
    # TODO: audit each line carefully
    chain = Dict{Block, Chain}()
    for j in transitives
        if j.left.block ∈ keys(chain) && j.right.block ∈ keys(chain)
            c₁, c₂ = chain[j.left.block], chain[j.right.block]
            if c₁ == c₂
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

Base.show(io::IO, G::Graph) = Base.show(io, (paths=values(G.sequence), blocks=values(G.block)))

# helper functions w/ common functionality
function write_fasta(io, name, seq)
    write(io, '>', name, '\n')
    write(io, columns(seq), '\n')
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

function marshal_fasta(io, G::Graph)
    for (i, b) in enumerate(values(G.block))
        write_fasta(io, b.uuid, b.sequence)
    end
end

# XXX: think of a way to break up function but maintain graph-wide node lookup table
function marshal_json(io, G::Graph)
    NodeID = NamedTuple{(:id,:name,:number,:strand), Tuple{String,String,Int,Bool}}
    nodes  = Dict{Node{Block}, NodeID}()

    # path serialization
    function dict(p::Path)
        blocks = Array{NodeID}(undef, length(p.node))
        counts = Dict{Block,Int}()

        for (i,node) ∈ enumerate(p.node)
            if node.block ∉ keys(counts)
                counts[node.block] = 1
            end
            blocks[i] = (
                id     = node.block.uuid,
                name   = p.name,
                number = counts[node.block], 
                strand = node.strand, 
            )
            nodes[node] = blocks[i]
            counts[node.block] += 1
        end

        return (
            name     = p.name,
            offset   = p.offset,
            circular = p.circular,
            blocks   = blocks,
        )
    end

    # block serialization
    pack(d::SNPMap) = [(k,Char(v)) for (k,v) ∈ d]
    pack(d::InsMap) = [(k,String(copy(v))) for (k,v) ∈ d]
    pack(d::DelMap) = [(k,v) for (k,v) ∈ d]

    function dict(b::Block)
        return (
            id       = b.uuid,
            sequence = String(sequence(b)),
            gaps     = b.gaps,
            mutate   = [(nodes[key], pack(val)) for (key,val) ∈ b.mutate],
            insert   = [(nodes[key], pack(val)) for (key,val) ∈ b.insert],
            delete   = [(nodes[key], pack(val)) for (key,val) ∈ b.delete]
        )
    end

    # NOTE: paths must come first as it fills the node lookup table
    paths  = [ dict(path)  for path  ∈ values(G.sequence) ]
    blocks = [ dict(block) for block ∈ values(G.block) ]

    JSON.print(io, (
        paths  = paths,
        blocks = blocks,
    ))
end

function marshal(io, G::Graph, fmt=:fasta)
    @match fmt begin
        :fasta || :fa => marshal_fasta(io, G)
        :json         => marshal_json(io, G)
        _ => error("$format not a recognized output format")
    end
end

# NOTE: only recognizes json input right now
function unmarshal(io)
    graph = JSON.parse(io)

    unpack = (
        snp = Dict(),
        ins = Dict(),
        del = Dict(),
    )
    blocks = Dict(map(graph["blocks"]) do blk
        # type wrangling
        b = (
            id       = String(blk["id"]),
            sequence = Array{UInt8}(blk["sequence"]),
            gaps     = Dict{Int,Int}(blk["gaps"]),
            mutate   = Dict(k=>v for (k,v) ∈ blk["mutate"]),
            insert   = Dict(k=>v for (k,v) ∈ blk["insert"]),
            delete   = Dict(k=>v for (k,v) ∈ blk["delete"]),
        )

        unpack.snp[b.id] = b.mutate
        unpack.ins[b.id] = b.insert
        unpack.del[b.id] = b.delete

        b.id => Block(
            b.id,
            b.sequence,
            b.gaps,
            # empty until we build the required node{block} objects
            Dict{Node{Block},SNPMap}(), 
            Dict{Node{Block},InsMap}(),
            Dict{Node{Block},DelMap}(),
        )
    end)

    paths = Dict(map(graph["paths"]) do path
        # type wrangling
        p = (
            name     = String(path["name"]),
            offset   = path["offset"],
            circular = path["circular"],
            blocks   = path["blocks"]
        )

        nodes = Node{Block}[]
        sizehint!(nodes, length(p.blocks))

        for (i,blk) ∈ enumerate(p.blocks)
            b = (
                id     = String(blk["id"]),
                strand = blk["strand"],
            )
            push!(nodes, Node(blocks[b.id], b.strand))

            # fill in block variant dictionaries
            blocks[b.id].mutate[nodes[i]] = Dict(
                snp[1] => UInt8(snp[2][1]) for snp ∈ unpack.snp[b.id][blk]
            )

            blocks[b.id].insert[nodes[i]] = Dict(
                ins[1] => Array{UInt8}(ins[2]) for ins ∈ unpack.ins[b.id][blk]
            )

            blocks[b.id].delete[nodes[i]] = Dict(
                del[1] => del[2] for del ∈ unpack.del[b.id][blk]
            )
        end
        p.name => Path(
            p.name,
            nodes,
            p.offset,
            p.circular,
        )
    end)

    return Graph(blocks, paths)
end

function test()
    graph = GZip.open("data/generated/assemblies/isolates.fna.gz", "r") do io
        isolates = graphs(io)
        println(">aligning...")
        align(isolates[1], isolates[2])
    end
    graph
end

end
