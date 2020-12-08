module Align

using FStrings, Match, Dates
using LinearAlgebra
using Infiltrator

import Base.Threads.@spawn

using ..Utility: read_paf, enforce_cutoff!
using ..Blocks
using ..Graphs

# Pool belongs only to the alignment module!
include("pool.jl")
using .Pool

export align

# ------------------------------------------------------------------------
# global variables
# TODO: move to a better location

fifos = pool(2)
getio()     = take!(fifos)
hasio()     = isready(fifos)
putio(fifo) = put!(fifos, fifo)
finalize()  = shutdown(fifos)
atexit(finalize)

# TODO: generalize to n > 2
function getios()
    @label getlock
    io₁ = getio()
    lock(fifos)
    if !hasio()
        putio(io₁)
        unlock(fifos)
        @goto getlock
    end
    io₂ = getio()
    unlock(fifos)

    return io₁, io₂
end

# ------------------------------------------------------------------------
# helper functions

function log(msg...)
    println(now(), " ", join(msg, " ")...)
    flush(stdout)
end

# command line execution
function execute(cmd::Cmd; now=true)
    out = Pipe()
    err = Pipe()

    proc = run(pipeline(ignorestatus(cmd), stdout=out, stderr=err); wait=now)

    close(out.in)
    close(err.in)

    stdout = @async String(read(out))
    stderr = @async String(read(err))

    if now
        return (
            out  = fetch(stdout),
            err  = fetch(stderr),
            code = proc.exitcode, #err,
        )
    else
        return (
            out  = stdout,
            err  = stderr,
            proc = proc,
        )
    end
end

function mash(input)
    result = execute(`mash triangle $input`)
    stdout = IOBuffer(result.out)

    N     = parse(Int64,readline(stdout))
    dist  = zeros(N,N)
    names = Array{String}(undef, N)
    for (i, line) in enumerate(eachline(stdout))
        elt = split(strip(line))
        names[i] = elt[1]
        dist[i,1:(i-1)] = [parse(Float64,x) for x in elt[2:end]]
    end

    dist = dist + dist';

    return dist, names
end

function minimap2(qry::String, ref::String)
    return execute(`minimap2 -x asm10 -m 10 -n 2 -s 30 -D -c $ref $qry`; now=false)
end

# ------------------------------------------------------------------------
# guide tree for order of pairwise comparison for multiple alignments

# TODO: distance?
mutable struct Clade
    name  ::String
    parent::Union{Clade,Nothing}
    left  ::Union{Clade,Nothing}
    right ::Union{Clade,Nothing}
    graph ::Channel{Graph}
end

# ---------------------------
# constructors

Clade()     = Clade("",nothing,nothing,nothing,Channel{Graph}(1))
Clade(name) = Clade(name,nothing,nothing,nothing,Channel{Graph}(1))
Clade(left::Clade, right::Clade) = Clade("",nothing,left,right,Channel{Graph}(1))

function Clade(distance, names; algo=:nj)
    @match algo begin
        :nj => return nj(distance, names)
        _ => error(f"unrecognized algorithm {algo}")
    end
end

function nj(distance, names)
    nodes = [Clade(names[i]) for i in 1:size(distance,1)]

    # internally scoped functions
    Q = function(D)
        n = size(D,1)
        q = (n-2)*D .- sum(D,dims=1) .- sum(D,dims=2)
        q[diagind(q)] .= Inf
        return q
    end

    pair = function(Q)
        ι  = argmin(Q)
        q₀ = Q[ι] 
        if ι.I[1] > ι.I[2]
            return ι.I[2], ι.I[1], q₀
        end
        return ι.I[1], ι.I[2], q₀
    end

    dist = function(D, i, j)
        n  = size(D,1)
        d₁ = .5*D[i,j] + 1/(2*(n-2))*(sum(D[i,:]) - sum(D[j,:])) 
        d₂ = D[i,j] - d₁

        if d₁ < 0
            d₂ -= d₁
            d₁ = 0
        end

        if d₂ < 0
            d₁ -= d₂
            d₂ = 0
        end

        dₙ = .5*(D[i,:] .+ D[j,:] .+ D[i,j])

        return d₁, d₂, dₙ
    end

    join! = function(D)
        q = Q(D)
        i, j, q₀ = pair(q)

        node = Clade(nodes[i], nodes[j])
        nodes[i].parent = node
        nodes[j].parent = node

        d₁, d₂, dₙ = dist(D, i, j)
        D[i,:] .= dₙ
        D[:,i] .= dₙ
        D[i,i]  = 0

        nodes[i] = node

        D = D[1:end .!= j, 1:end .!= j]
        deleteat!(nodes, j)

        return D
    end

    # body of function
    @assert length(names) == length(Set(names))

    while length(nodes) > 2
        distance = join!(distance)
    end

    return Clade(nodes[1], nodes[2])
end

# ---------------------------
# operators

isleaf(c::Clade) = isnothing(c.left) && isnothing(c.right)

# serialization to newick format
function Base.show(io::IO, c::Clade) 
    if isleaf(c)
        print(io, c.name)
    else
        print(io, "(")
        show(io, c.left)
        print(io, ",")
        show(io, c.right)
        print(io, ")")
    end
end

function leaves(root::Clade)
    itr = Channel{Clade}(0)
    function traverse(node::Clade)
        if isleaf(node)
            put!(itr, node)
        else
            traverse(node.left)
            traverse(node.right)
        end
    end

    @async begin
        traverse(root)
        close(itr)
    end

    return itr
end

function postorder(root::Clade)
    itr = Channel{Clade}(0)
    function traverse(node::Clade)
        if isleaf(node)
            put!(itr, node)
        else
            traverse(node.left)
            traverse(node.right)
            put!(itr,node)
        end
    end

    @async begin
        traverse(root)
        close(itr)
    end

    return itr
end

function preorder(root::Clade)
    itr = Channel{Clade}(0)
    function traverse(node::Clade)
        if isleaf(node)
            put!(itr, node)
        else
            put!(itr,node)
            traverse(node.left)
            traverse(node.right)
        end
    end

    @async begin
        traverse(root)
        close(itr)
    end

    return itr
end

# ------------------------------------------------------------------------
# ordering functions

# TODO: assumes the input graphs are singletons! generalize
function ordering(Gs...; compare=mash)
    fifo = getio()

    task = @async compare(path(fifo))

    @spawn begin
        open(fifo,"w") do io
            for G in Gs
                serialize(io, G)
            end
        end
    end

    distance, names = fetch(task)
    root = Clade(distance, names; algo=:nj)

    putio(fifo)
    return root
end

# ------------------------------------------------------------------------
# align functions

function align_pair(G₁::Graph, G₂::Graph, energy::Function)
    function write(fifo, G)
        open(fifo, "w") do io
            marshal(io, G)
        end
    end

    io₁, io₂ = getios()

    cmd = minimap2(path(io₁), path(io₂))

    # NOTE: minimap2 opens up file descriptors in order!
    #       must process 2 before 1 otherwise we deadlock
    write(io₂, G₂) # ref
    write(io₁, G₁) # qry

    out  = IOBuffer(fetch(cmd.out))
    hits = collect(read_paf(out))
    sort!(hits; by=energy)

    close(out)
    putio(io₁)
    putio(io₂)

    # NOTE: we could turn this section into its own function
    blocks = Dict{String,Block}()
    for hit in hits
        log(hit)
        if energy(hit) >= 0
            break
        end

        if !(hit.qry.name in keys(G₁.block)) || !(hit.ref.name in keys(G₂.block))
            continue
        end

        qry₀ = pop!(G₁.block, hit.qry.name)
        ref₀ = pop!(G₂.block, hit.ref.name)

        hit.qry.seq = qry₀.sequence
        hit.ref.seq = ref₀.sequence

        enforce_cutoff!(hit, 100) # TODO: remove hard-coded parameter

        blks = combine(qry, ref, hit)
        qrys = filter(b -> b.kind != :ref, blks)
        refs = filter(b -> b.kind != :qry, blks)

        for (isolate, path) in G₁.sequence
            replace!(path, qry₀, qrys)
        end

        for (isolate, path) in G₂.sequence
            replace!(path, ref₀, refs)
        end

        for blk in blks
            blocks[blk.uuid] = blk
        end
    end

    sequence = merge(G₁.sequence, G₂.sequence)

    # XXX: worry about uuid collision?
    merge!(blocks, G₁.block)
    merge!(blocks, G₂.block)

    # TODO: remove transitives
    
    return Graph(
            blocks, 
            sequence,
            # ...
    )
end

# TODO: the associate array is a bit hacky...
#       can we push it directly into the channel?
function align(Gs::Graph...; energy=(hit)->(-Inf))
    function kernel(clade)
        Gₗ = take!(clade.left.graph)
        Gᵣ = take!(clade.right.graph)
        G₀ = align_pair(Gₗ, Gᵣ, energy)

        # TODO: self-maps!

        put!(clade.graph, G₀)
    end

    log("--> ordering")
    tree = ordering(Gs...)
    tips = Dict{String,Graph}(collect(keys(G.sequence))[1] => G for G in Gs)

    for clade in postorder(tree)
        if isleaf(clade)
            put!(clade.graph, tips[clade.name])
            close(clade.graph)
        else
            @sync begin
                kernel(clade)
                close(clade.graph)
            end
        end
    end

    return take!(tree.graph)
end

# ------------------------------------------------------------------------
# testing

function test()
    log("> testing mash...")
    distance, names = mash("data/generated/assemblies/isolates.fna.gz")
    tree = Clade(distance, names; algo=:nj) 
    log("done!")
end

end
