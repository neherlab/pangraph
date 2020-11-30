module Align

using FStrings, Match
using LinearAlgebra

# NOTE: this is temporary
include("pool.jl")
using .Pool

include("graph.jl")
using .Pangraph

# ------------------------------------------------------------------------
# global variables
# TODO: move to a better location

fifos = pool(10)
getio()     = take!(fifos)
putio(fifo) = put!(fifos, fifo)

# ------------------------------------------------------------------------
# helper functions

# command line execution
function execute(cmd::Cmd)
    out = Pipe()
    err = Pipe()

    proc = run(pipeline(ignorestatus(cmd), stdout=out, stderr=err))

    close(out.in)
    close(err.in)

    stdout = @async String(read(out))
    stderr = @async String(read(err))

    wait(proc)
    return (
        out  = fetch(stdout),
        err  = fetch(stderr),
        code = proc.exitcode
    )
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
    return execute(`minimap2 -x asm20 -m 10 -n 2 -s 30 -D -c $ref $qry`)
end

# ------------------------------------------------------------------------
# guide tree for order of pairwise comparison for multiple alignments

# TODO: distance?
mutable struct Clade
    name  ::String
    parent::Union{Clade,Nothing}
    left  ::Union{Clade,Nothing}
    right ::Union{Clade,Nothing}
end

# ---------------------------
# constructors

Clade()     = Clade("",nothing,nothing,nothing)
Clade(name) = Clade(name,nothing,nothing,nothing)
Clade(left::Clade, right::Clade) = Clade("",nothing,left,right)

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

    return Clade(nodes[1], nodes[2]) end

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
function ordering(Gs...)
    fifo = getio()
    @spawn begin
        open(path(fifo)) do io
            for G in Gs
                write(io, G)
            end
        end
    end

    # TODO: allow for other pairwise comparisons besides mash
    distance, names = mash(path(fifo))
    root = Clade(distance, names; algo=:nj)

    putio(fifo)

    return root
end

# ------------------------------------------------------------------------
# align functions

module Minimap2

function merge(G1::Graph, G2::Graph)
    return G1
end

# TODO: the associate array is a bit hacky...
function align(Gs::Graph...)
    tips  = Dict([(keys(G.sequence)[1],G) for G in Gs])
    tree  = ordering(Gs...)

    merge = Dict{Clade,Graph}()
    for clade in tree
        if isleaf(clade)
            merged[clade] = tips[clade.name]
            continue
        end

        merged[clade] = merge(merged[clade.left], merged[clade.right])

        delete!(merged, clade.left)
        delete!(merged, clade.right)
    end
end

end

# ------------------------------------------------------------------------
# testing

function test()
    println("testing mash...")
    distance, names = mash("data/generated/assemblies/isolates.fna.gz")
    tree = Clade(distance, names; algo=:nj) 
    println(tree)
    println("done!")

    println("testing alignment...")

    shutdown(fifos)
end

end
