module Align

using FStrings, Match, Dates
using LinearAlgebra
using Infiltrator

import GZip # NOTE: only for debugging
import Base.Threads.@spawn

# NOTE: this is temporary
include("util.jl")
using .Utility: read_paf

include("pool.jl")
using .Pool

include("graph.jl")
using .Pangraph

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
function getios(i)
    @label getlock
    log("-----> getting io 1", i)
    io₁ = getio()
    lock(fifos)
    if !hasio()
        putio(io₁)
        log("-----> resetting io 1", i)
        unlock(fifos)
        @goto getlock
    end
    log("-----> getting io 2", i)
    io₂ = getio()
    log("-----> obtained ios", i, path(io₁), path(io₂))
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

    log("-------> pipeline...")
    proc = run(pipeline(ignorestatus(cmd), stdout=out, stderr=err); wait=now)

    log("-------> closing...")
    close(out.in)
    close(err.in)

    log("-------> output...")
    stdout = @async String(read(out))
    stderr = @async String(read(err))

    if now
        log("-------> fetching...")
        return (
            out  = fetch(stdout),
            err  = fetch(stderr),
            code = proc.exitcode, #err,
        )
    else
        log("-------> started...", process_running(proc))
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
    return execute(`minimap2 -x asm20 -m 10 -n 2 -s 30 -D -c $ref $qry`; now=false)
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

# TODO: need to fix. skips right branch of root
function Base.iterate(root::Clade, state)
    node, last = state

    if isnothing(node)
        return nothing
    end

    if isleaf(node)
        @goto yield
    end

    if isnothing(last) || last == root
        return iterate(node, (node.left, node))
    end

    if last == node.left
        return iterate(node, (node.right, node))
    end

    if last == node.right
        @goto yield
    end

    # if we reach here, we have a serious problem!
    error("invalid iteration case")

    @label yield
    return node, (node.parent, node)
end

function Base.iterate(root::Clade)
    return iterate(root, (root, nothing))
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

# TODO: fill in actual implementation
function merge(G₁::Graph, G₂::Graph, i)
    function write(fifo, G)
        open(fifo, "w") do io
            marshal(io, G)
        end
    end
    log("-----> grabbing fifos", i)
    io₁, io₂ = getios(i)

    log("-----> starting minimap2", i)
    cmd = minimap2(path(io₁), path(io₂))

    log("-----> opening ios", i)
    write(io₂, G₂)
    log("-----> wrote 2", i)
    write(io₁, G₁)

    log("-----> parsing minimap2", i)
    hits = read_paf(IOBuffer(fetch(cmd.out)))

    log("-----> putting io 1", i)
    putio(io₁)
    log("-----> putting io 2", i)
    putio(io₂)

    return G₁
end

# TODO: the associate array is a bit hacky...
#       can we push it directly into the channel?
function align(Gs::Graph...)
    function kernel(i, clade)
        log("--> waiting for children")
        Gₗ = take!(clade.left.graph)
        Gᵣ = take!(clade.right.graph)
        log("--> merging", i, clade == tree, clade == tree.left, clade == tree.right)
        put!(clade.graph, merge(Gₗ, Gᵣ, i))
    end

    tree = ordering(Gs...)
    tips = Dict{String,Graph}(collect(keys(G.sequence))[1] => G for G in Gs)

    log("--> aligning...")
    for (i, clade) in enumerate(postorder(tree))
        if isleaf(clade)
            put!(clade.graph, tips[clade.name])
        else
            kernel(i, clade)
        end
        close(clade.graph)
    end

    log("--> waiting on root")
    return take!(tree.graph)
end

# ------------------------------------------------------------------------
# testing

function test()
    log("> testing mash...")
    distance, names = mash("data/generated/assemblies/isolates.fna.gz")
    tree = Clade(distance, names; algo=:nj) 
    log("done!")

    log("> testing alignment...")
    GZip.open("data/generated/assemblies/isolates.fna.gz", "r") do io
        graphs = Graphs(io)
        graph  = align(graphs...)
    end
end

end
