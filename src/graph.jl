module Graphs

using GZip # NOTE: for debugging purposes
using Random
using Rematch

import JSON

# ---------------------------
# functions to extend in submodules

export pair
function pair(item) end

export sequence
function sequence(obj, name; gaps=false)     end
function sequence(obj; gaps=false)           end
function sequence!(s, obj, name; gaps=false) end
function sequence!(s, obj; gaps=false)       end

export reverse_complement
function reverse_complement(item) end

include("interval.jl")
include("counter.jl")
include("util.jl")
include("node.jl")
include("block.jl")
include("path.jl")
include("junction.jl")

using .Utility: read_fasta, name, columns, log
using .Nodes
using .Blocks
using .Paths
using .Junctions
using .Intervals

export Graph
export graphs, marshal, serialize, detransitive!

# ------------------------------------------------------------------------
# graph data structure

struct Graph
    block    :: Dict{String,Block}   # uuid      -> block
    sequence :: Dict{String,Path}    # isolation -> path
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

graphs(io::IO; circular=false) = [Graph(record.name, record.seq; circular=circular) for record in read_fasta(io)]

# --------------------------------
# operators

const Link  = NamedTuple{(:block,:strand),Tuple{Block, Bool}}
const Chain = Array{Link, 1}

# XXX: break into smaller functions
#      too long
function detransitive!(G::Graph)
    isosᵥ = count_isolates(values(G.sequence))

    # collect all junctions that transitively pass isolates through
    transitives = Junction[]
    for (j, isosⱼ) in junctions(values(G.sequence))
        if ((isosᵥ[j.left.block] == isosᵥ[j.right.block]) &&
            (isosᵥ[j.left.block] == isosⱼ))
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

            c₁ == c₂ && continue

            merged =
            if left(j) == last(c₁) && right(j) == first(c₂)
                cat(c₁, c₂, dims=1)
            elseif left(j) == last(c₁) && rev(right(j)) == last(c₂)
                cat(c₁, rev(c₂), dims=1)
            elseif rev(left(j)) == first(c₁) && right(j) == first(c₂)
                cat(rev(c₁), c₂, dims=1)
            elseif rev(left(j)) == first(c₁) && rev(right(j)) == last(c₂)
                cat(c₂, c₁, dims=1)
            else
                error("case not covered")
            end

            for b in first.(merged)
                chain[b] = merged
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
        @assert all([isosᵥ[C.block] == isos for C in c[2:end]])
        new  = Block((s ? b : reverse_complement(b) for (b,s) ∈ c)...)

        for iso ∈ keys(isos)
            replace!(G.sequence[iso], c, new)
        end

        for b ∈ first.(c)
            pop!(G.block, b.uuid)
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

# ------------------------------------------------------------------------
# operators

function sequence(g::Graph, name::AbstractString)
    name ∉ keys(g.sequence) && error("'$name' not a valid sequence identifier")
    path = g.sequence[name]
    return join(sequence(node.block, node) for node ∈ path.node)
end

sequence(g::Graph) = [ name => join(sequence(node.block, node) for node ∈ path.node) for (name, path) ∈ g.sequence ]

# ------------------------------------------------------------------------
# main point of entry

using Random: seed!

function test(file="data/generated/assemblies/isolates.fna.gz")
    seed!(0)

    open = endswith(file,".gz") ? GZip.open : Base.open

    log("> running block test...")
    if !Blocks.test()
        error("failed individual block reconstruction")
    end

    index = 1:50
    log("> running graph test...")
    log("-> building graph...")
    graph, isolates = open(file, "r") do io
        isolates = graphs(io)
        println(">aligning...")
        align(isolates[index]...;maxgap=10) , isolates
    end

    log("-> verifying graph...")
    for isolate ∈ isolates[index]
        name, seq₀ = first(sequence(isolate))
        seq₁ = sequence(graph, name)
        if !all(seq₀ .== seq₁)
            error("incorrect sequence reconstruction")
        else
            log("---> isolate '$name' correctly reconstructed")
        end
    end

    graph
end

end
