module Graphs

using GZip # NOTE: for debugging purposes
using Random
using Rematch
using ProgressMeter

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

export reverse_complement, reverse_complement!
function reverse_complement(item)  end
function reverse_complement!(item) end

export marshal, marshal_fasta, marshal_json
function marshal_fasta(io::IO, x) end
function marshal_json(io::IO, x) end
function marshal(io::IO, x, fmt=:fasta)
    @match fmt begin
        :fasta || :fa => return marshal_fasta(io, x)
        :json         => return marshal_json(io, x)
        _ => error("$fmt not a recognized output format")
    end
end

export serialize
function serialize(io::IO, x) end

# ------------------------------------------------------------------------
# aux types

Maybe{T} = Union{T,Nothing}

# aliases
const SNPMap = Dict{Int,UInt8}
const InsMap = Dict{Tuple{Int,Int},Array{UInt8,1}}
const DelMap = Dict{Int,Int}

export Maybe, SNPMap, InsMap, DelMap

Base.show(io::IO, m::SNPMap) = show(io, [ k => Char(v) for (k,v) in m ])
Base.show(io::IO, m::InsMap) = show(io, [ k => String(Base.copy(v)) for (k,v) in m ])

include("interval.jl")
include("counter.jl")
include("node.jl")
include("util.jl")
include("block.jl")
include("path.jl")
include("edge.jl")
# include("pool.jl")
include("junction.jl")
include("cmd.jl")
include("minimap.jl")
include("mash.jl")

using .Utility:
    read_fasta, write_fasta, name, columns, log,
    make_consensus, alignment_alleles
using .Nodes
using .Blocks
using .Paths
using .Edges
using .Junctions
using .Intervals
# using .Pool

import .Shell: mash, mafft, havecommand
import .Minimap: PanContigs

export Graph

export graphs, detransitive!, prune!, finalize!
export pancontigs
export checkblocks

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
    numisos = count_isolates(values(G.sequence))

    # collect all junctions that transitively pass isolates through
    transitives = Junction[]
    for (j, depth) in junctions(values(G.sequence))
        if (numisos[j.left.block] == numisos[j.right.block] == depth)
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
        isos = numisos[c[1].block]
        @assert all([numisos[C.block] == isos for C in c[2:end]])

        new = Block((s ? b : reverse_complement(b) for (b,s) ∈ c)...)
        for iso ∈ keys(isos)
            # oldseq = sequence(G.sequence[iso])

            replace!(G.sequence[iso], c, new)

            # newseq = sequence(G.sequence[iso])
            #= if oldseq != newseq
                path = G.sequence[iso]
                badloci = Int[]
                for i ∈ 1:min(length(newseq),length(oldseq))
                    if newseq[i] != oldseq[i]
                        push!(badloci, i)
                    end
                end
                left, right = max(badloci[1]-10, 1), min(badloci[1]+10, length(newseq))
                cumulative_lengths = cumsum([length(n.block, n) for n in path.node])

                println("--> length:           ref($(length(oldseq))) <=> seq($(length(newseq)))")
                println("--> number of nodes:  $(length(path.node))")
                println("--> |badloci|:        $(length(badloci))")
                println("--> window:           $(left):$(badloci[1]):$(right)")
                println("--> ref:              $(oldseq[left:right])")
                println("--> seq:              $(newseq[left:right])")
            end =#
        end

        for b ∈ first.(c)
            pop!(G.block, b.uuid)
        end

        G.block[new.uuid] = new

        #=
        for b ∈ values(G.block)
            check(b)
        end
        =#
    end
    # checkblocks(G)
end

function prune!(G::Graph)
    used = Set(n.block.uuid for p in values(G.sequence) for n in p.node)
    filter!((blk)->first(blk) ∈ used, G.block)
end

function checkblocks(G::Graph)
    used   = Set(n.block.uuid for p in values(G.sequence) for n in p.node)
    stored = Set(keys(G.block))

    if used != stored
        @show setdiff(used,stored)
        @show setdiff(stored,used)
        # @infiltrate
        error("bad blocks")
    end
end

# ------------------------------------------------------------------------
# i/o & (de)serialization

Base.show(io::IO, G::Graph) = Base.show(io, (paths=values(G.sequence), blocks=values(G.block)))

pancontigs(G::Graph) = let
    uuid = collect(Base.keys(G.block))
    PanContigs(
        uuid,
        [String(sequence(G.block[id])) for id in uuid],
    )
end

# TODO: can we generalize to multiple individuals
#       equivalent to "highway" detection
function serialize(io::IO, G::Graph)
    if length(G.sequence) != 1
        error("only singleton graphs implemented")
    end

    name = collect(keys(G.sequence))[1]
    seq  = collect(values(G.block))[1].sequence

    write_fasta(io, name, seq)
end

function marshal_fasta(io::IO, G::Graph)
    for (i, b) in enumerate(values(G.block))
        write_fasta(io, b.uuid, b.sequence)
    end
end

# XXX: think of a way to break up function but maintain graph-wide node lookup table
function marshal_json(io::IO, G::Graph)
    NodeID    = NamedTuple{(:id,:name,:number,:strand), Tuple{String,String,Int,Bool}}
    nodes     = Dict{Node{Block}, NodeID}()
    positions = Dict{Block, Dict{NodeID, Tuple{Int,Int}}}()

    finalize!(G)

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

            if node.block ∉ keys(positions)
                positions[node.block] = Dict{NodeID, Tuple{Int,Int}}()
            end
            positions[node.block][blocks[i]] = (p.position[i], p.position[i == length(p.node) ? 1 : i+1] - 1)
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

    strip(id) = (name=id.name,number=id.number,strand=id.strand)
    function dict(b::Block)
        return (
            id        = b.uuid,
            sequence  = String(sequence(b)),
            gaps      = b.gaps,
            mutate    = [(strip(nodes[key]), pack(val)) for (key,val) ∈ b.mutate],
            insert    = [(strip(nodes[key]), pack(val)) for (key,val) ∈ b.insert],
            delete    = [(strip(nodes[key]), pack(val)) for (key,val) ∈ b.delete],
            positions = [(strip(key), val) for (key,val) ∈ positions[b]]
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
            gaps     = Dict{Int,Int}(
                parse(Int,k) => v for (k,v) ∈ blk["gaps"]
            ),
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
            filter!(p->p.first!="id", blk)

            blocks[b.id].mutate[nodes[i]] = Dict(
                snp[1] => UInt8(snp[2][1]) for snp ∈ unpack.snp[b.id][blk]
            )

            blocks[b.id].insert[nodes[i]] = Dict(
                let
                    @assert length(ins[1]) == 2
                    key = (ins[1][1], ins[1][2])
                    key => Array{UInt8}(ins[2])
                end for ins ∈ unpack.ins[b.id][blk]
            )

            blocks[b.id].delete[nodes[i]] = Dict(
                del[1] => del[2] for del ∈ unpack.del[b.id][blk]
            )
        end

        path = Path(
            p.name,
            nodes,
            p.offset,
            p.circular,
            [],
        )
        positions!(path)

        p.name => path
    end)

    return Graph(blocks, paths)
end

# ------------------------------------------------------------------------
# operators

function sequence(g::Graph, name::AbstractString)
    name ∉ keys(g.sequence) && error("'$name' not a valid sequence identifier")
    path = g.sequence[name]

    return sequence(path)
end

sequence(g::Graph) = [ name => join(String(sequence(node.block, node)) for node ∈ path.node) for (name, path) ∈ g.sequence ]

function realign!(g::Graph)
    meter = Progress(length(g.block); desc="polishing progress", output=stderr)
    Threads.@threads for blk in collect(values(g.block))
        io, node = mafft(blk)

        seq = collect(read_fasta(io))
        aln = reduce(hcat, map((r)->r.seq, seq))
        ref = make_consensus(aln)

        nodes = map((r)->node[r.name], seq)
        blk.gaps, blk.mutate, blk.delete, blk.insert, blk.sequence = alignment_alleles(ref, aln, nodes)
        next!(meter)
    end
end

function finalize!(g)
    for p in values(g.sequence)
        positions!(p)
    end

    # realign!(g)
end

# ------------------------------------------------------------------------
# main point of entry

function test(file="data/marco/mycobacterium_tuberculosis/genomes.fa")
    open = endswith(file,".gz") ? GZip.open : Base.open

    log("> running graph test...")
    log("-> building graph...")

    sequences = String[]
    energy = (aln) -> let
        len = aln.length
        len < 100 && return Inf

        cuts(hit) = (hit.start > 100) + ((hit.length-hit.stop) > 100)

        ncuts = cuts(aln.qry)+cuts(aln.ref)
        nmuts = aln.divergence*aln.length

        return -len + 100*ncuts + 20*nmuts
    end

    graph, isolates = open(file, "r") do io
        isolates  = graphs(io; circular=true)
        sequences = [first(sequence(iso)) for iso in isolates]

        println("-->aligning...")
        align(isolates...;energy=energy,minblock=100), isolates
    end

    log("-> verifying graph...")
    for (i, isolate) ∈ enumerate(isolates)
        name, seq₀ = sequences[i]
        seq₁ = sequence(graph, name)
        if !all(seq₀ .== seq₁)
            path = graph.sequence[name]
            x    = findfirst(seq₁, "ACTTGGCTATCCCGCAGGAC")

            println("> true ($(length(seq₀))):          ", seq₀[1:20])
            println("> reconstructed ($(length(seq₁))): ", seq₁[1:20])
            println("> offset:                          $(path.offset)")
            println("> needed offset:                   $(x)")
            log("--> isolate '$name' incorrectly reconstructed")
        else
            log("--> isolate '$name' correctly reconstructed")
        end
    end

    finalize!(graph)

    graph
end

end
