module Align

using Rematch, Dates
using LinearAlgebra
using ProgressMeter

using Base.Threads: @spawn, @threads

using ..Utility: read_paf, enforce_cutoff!
using ..Blocks
using ..Paths: replace!
using ..Nodes
using ..Graphs
using ..Mash

import ..Minimap

export align

# ------------------------------------------------------------------------
# helper functions

function log(msg...)
    println(stderr, now(), " ", join(msg, " ")...)
    flush(stderr)
end

# ------------------------------------------------------------------------
# guide tree for order of pairwise comparison for multiple alignments

# TODO: distance?
"""
	mutable struct Clade
		name   :: String
		parent :: Union{Clade,Nothing}
		left   :: Union{Clade,Nothing}
		right  :: Union{Clade,Nothing}
		graph  :: Channel{Graph}
	end

Clade is a node (internal or leaf) of a binary guide tree used to order pairwise alignments
associated to a multiple genome alignment in progress.
`name` is only non-empty for leaf nodes.
`parent` is `nothing` for the root node.
`graph` is a 0-sized channel that is used as a message passing primitive in alignment.
"""
mutable struct Clade
    name   :: String
    parent :: Union{Clade,Nothing}
    left   :: Union{Clade,Nothing}
    right  :: Union{Clade,Nothing}
    graph  :: Channel{Graph}
end

# ---------------------------
# constructors

"""
	Clade()

Generate an empty, disconnected clade.
"""
Clade()     = Clade("",nothing,nothing,nothing,Channel{Graph}(2))
"""
	Clade(name)

Generate an empty, disconnected clade with name `name`.
"""
Clade(name) = Clade(name,nothing,nothing,nothing,Channel{Graph}(2))
"""
	Clade(left::Clade, right::Clade)

Generate an nameless clade with `left` and `right` children.
"""
function Clade(left::Clade, right::Clade)
    parent = Clade("",nothing,left,right,Channel{Graph}(2))

    left.parent = parent
    right.parent = parent

    parent
end

"""
	Clade(distance, names; algo=:nj)

Generate a tree from a matrix of pairwise distances `distance`.
The names of leafs are given by an array of strings `names`.
`algo` dictates the algorithm used to transform the distance matrix into a tree.
Currently on neighbor joining (:nj) is supported.
"""
function Clade(distance, names; algo=:nj)
    @match algo begin
        :nj => return nj(distance, names)
        _ => error("unrecognized algorithm $(algo)")
    end
end

"""
	nj(distance, names)

Lower-level function.
Generate a tree from a matrix of pairwise distances `distance`.
The names of leafs are given by an array of strings `names`.
Uses neighbor joining.
"""
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

"""
	isleaf(c::Clade)

Return if Clade `c` is a terminal node, i.e. a leaf.
"""
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

"""
	leaves(root::Clade)

Return all terminal nodes that have `root` as their common ancestor.
"""
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

"""
	postorder(root::Clade)

Return an postorder iterator over descendents of `root`.
"""
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

"""
	preorder(root::Clade)

Return an preorder iterator over descendents of `root`.
"""
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

function bisect(tips)
    length(tips) > 1 || return tips[1]

    m = length(tips) ÷ 2
    l = bisect(tips[1:m])
    r = bisect(tips[(m+1):end])

    return Clade(l, r)
end

"""
	balance(root::Clade)

Rotate a binary tree to a balanced configuration.
Preserves the topological ordering of the leaves.
"""
function balance(root::Clade)
    tips = collect(leaves(root))
    return bisect(tips)
end

function Base.length(node::Clade)
    isleaf(node) && return 1
    return 1 + Base.length(node.left) + Base.length(node.right)
end

# ------------------------------------------------------------------------
# ordering functions

# TODO: assumes the input graphs are singletons! generalize
"""
	ordering(compare, Gs...)

Return a guide tree based upon distances computed from a collection of graphs `Gs`, using method `compare`.
The signature of `compare` is expected to be `compare(G::Graphs....) -> distance, names`.
Graphs `Gs...` are expected to be singleton graphs.
"""
function ordering(compare, Gs...)
    distance, names = compare(Gs...)
    return Clade(distance, names; algo=:nj)
end

# ------------------------------------------------------------------------
# align functions

function preprocess(hits, skip, energy, blocks!)
    hits = [
        let
            qry₀, ref₀ = blocks!(hit)

            hit.qry.seq = qry₀.sequence
            hit.ref.seq = ref₀.sequence

            qry, strand =
            if !hit.orientation
                reverse_complement!(hit)
                reverse_complement(qry₀; keepid=true), false
            else
                qry₀, true
            end
            ref = ref₀

            (
                hit=hit,
                qry=qry,
                ref=ref,
                qry₀=qry₀,
                ref₀=ref₀,
                strand=strand,
           )
        end for hit in hits if (energy(hit) < 0 && !skip(hit))
    ]

    return hits
end

function do_align(G₁::Graph, G₂::Graph, energy::Function, minblock::Int, sensitivity::String)
    hits = Minimap.align(pancontigs(G₁), pancontigs(G₂), minblock, sensitivity)
    sort!(hits; by=energy)

    return hits
end

# TODO: make thread safe
function align_kernel(matches, minblock, replace!, verbose)
    blocks   = Array{Array{Block,1},1}(undef, length(matches))
    for (i, match) in enumerate(matches)
        # destructure precomputed hit
        hit,qry,ref,qry₀,ref₀,strand = match

        verbose && log(hit)

        enforce_cutoff!(hit, minblock)
        blks = combine(qry, ref, hit; minblock=minblock)

        qrys = map(b -> b.block, filter(b -> b.kind != :ref, blks))
        refs = map(b -> b.block, filter(b -> b.kind != :qry, blks))

        replace!((qry=qry₀, ref=ref₀), (qry=qrys, ref=refs), strand)

        blocks[i] = map(b->b.block, blks)
    end

    return Dict(blk.uuid => blk for block in blocks for blk in block)
end

"""
	align_self(G₁::Graph, energy::Function, minblock::Int, verify::Function, verbose::Bool; maxiter=100)

Align graph `G₁` to itself by looking for homology between blocks.
This is a low-level function.

`energy` is to be a function that takes an alignment between two blocks and produces a score.
The _lower_ the score, the _better_ the alignment. Only negative energies are considered.

`minblock` is the minimum size block that will be produced from the algorithm.
`maxiter` is maximum number of duplications that will be considered during this alignment.
"""
function align_self(G₁::Graph, energy::Function, minblock::Int, verify::Function, verbose::Bool; maxiter=100, sensitivity="asm10")
    G₀ = G₁

    for niter in 1:maxiter
        # calculate pairwise hits
        hits = do_align(G₀, G₀, energy, minblock, sensitivity)

        skip = (hit) -> (
               (hit.qry.name == hit.ref.name)
            || (hit.length < minblock)
            || (!(hit.qry.name in keys(G₀.block)) || !(hit.ref.name in keys(G₀.block)))
        )
        block = (hit) -> (
            qry = pop!(G₀.block, hit.qry.name),
            ref = pop!(G₀.block, hit.ref.name),
        )
        replace = (old, new, orientation) -> let
            for path in values(G₀.sequence)
                replace!(path, old.qry, new.qry, orientation)
            end

            for path in values(G₀.sequence)
                replace!(path, old.ref, new.ref, true)
            end
        end

        # filter out according to energy & skip conditions
        # only continue if at least one hit survived
        merges = preprocess(hits, skip, energy, block)
        length(merges) > 0 || break

        blocks = align_kernel(merges, minblock, replace, verbose)
        merge!(blocks, G₀.block)

        G₀ = Graph(
            blocks,
            G₀.sequence,
        )
        detransitive!(G₀)
        purge!(G₀)
        prune!(G₀)
    end

    return G₀
end

function compare(old, new, strand)
    for node in keys(old)
        oldseq = sequence(old, node)

        δ         = 1
        newblkseq = UInt8[]
        oldblkseq = UInt8[]
        oldstrand = node.strand

        if oldstrand && strand
            for blk in new
                newblkseq = sequence(blk, node)
                oldblkseq = oldseq[δ:δ+length(newblkseq)-1]
                if !all(oldblkseq .== newblkseq)
                    @goto bad
                end

                δ += length(newblkseq)
            end
        elseif oldstrand && !strand
            for blk in reverse(new)
                newblkseq = reverse_complement(sequence(blk, node))
                oldblkseq = oldseq[δ:δ+length(newblkseq)-1]
                if !all(oldblkseq .== newblkseq)
                    @goto bad
                end

                δ += length(newblkseq)
            end
        elseif !oldstrand && strand
            for blk in reverse(new)
                newblkseq = sequence(blk, node)
                oldblkseq = oldseq[δ:δ+length(newblkseq)-1]
                if !all(oldblkseq .== newblkseq)
                    @goto bad
                end

                δ += length(newblkseq)
            end
        else
            for blk in new
                newblkseq = reverse_complement(sequence(blk, node))
                oldblkseq = oldseq[δ:δ+length(newblkseq)-1]
                if !all(oldblkseq .== newblkseq)
                    @goto bad
                end

                δ += length(newblkseq)
            end
        end

        continue

    @label bad
        badblkloci = Int[]
        for i ∈ 1:min(length(newblkseq),length(oldblkseq))
            if newblkseq[i] != oldblkseq[i]
                push!(badblkloci, i)
            end
        end
        left, right = max(badblkloci[1]-10, 1), min(badblkloci[1]+10, length(newblkseq))

        println("--> length:           ref($(length(oldblkseq))) <=> seq($(length(newblkseq)))")
        println("--> number bad:       $(length(badblkloci))")
        println("--> window:           $(left):$(badblkloci[1]):$(right)")
        println("--> old:              $(String(oldblkseq[left:right]))") 
        println("--> new:              $(String(newblkseq[left:right]))") 

        @show strand, oldstrand

        error("bad sequence reconstruction")
    end
end

"""
	align_pair(G₁::Graph, G₂::Graph, energy::Function, minblock::Int, verify::Function, verbose::Bool; maxiter=100)


Align graph `G₁` to graph `G₂` by looking for homology between consensus sequences of blocks.
This is a low-level function.

`energy` is to be a function that takes an alignment between two blocks and produces a score.
The _lower_ the score, the _better_ the alignment. Only negative energies are considered.

`minblock` is the minimum size block that will be produced from the algorithm.
`maxiter` is maximum number of duplications that will be considered during this alignment.
"""
function align_pair(G₁::Graph, G₂::Graph, energy::Function, minblock::Int, verify::Function, verbose::Bool; sensitivity="asm10")
    hits = do_align(G₁, G₂, energy, minblock, sensitivity)

    skip = (hit) -> (
            !(hit.qry.name in keys(G₁.block))
         || !(hit.ref.name in keys(G₂.block))
         || (hit.length < minblock)
    )
    block = (hit) -> (
        qry = pop!(G₁.block, hit.qry.name),
        ref = pop!(G₂.block, hit.ref.name),
    )

    replace = (old, new, orientation) -> let
        for path in values(G₁.sequence)
            replace!(path, old.qry, new.qry, orientation)
        end

        for path in values(G₂.sequence)
            replace!(path, old.ref, new.ref, true)
        end
    end

    merges = preprocess(hits, skip, energy, block)

    blocks   = align_kernel(merges, minblock, replace, verbose)
    sequence = merge(G₁.sequence, G₂.sequence)

    # XXX: worry about uuid collision?
    merge!(blocks, G₁.block)
    merge!(blocks, G₂.block)

    G = Graph(
        blocks,
        sequence,
    )

    detransitive!(G)
    purge!(G)
    prune!(G)

    return G
end

# TODO: the associative array is a bit hacky...
#       can we push it directly into the channel?
"""
	align(G::Graph...; compare=Mash.distance, energy=(hit)->(-Inf), minblock=100)

Align graph `G₁` to itself by looking for homology between blocks.
Multithreaded by default.
This is usually the function you want.

`energy` is to be a function that takes an alignment between two blocks and produces a score.
The _lower_ the score, the _better_ the alignment. Only negative energies are considered.

`minblock` is the minimum size block that will be produced from the algorithm.
`maxiter` is maximum number of duplications that will be considered during this alignment.

`compare` is the function to be used to generate pairwise distances that generate the internal guide tree.
"""
function align(Gs::Graph...; compare=Mash.distance, energy=(hit)->(-Inf), minblock=100, reference=nothing, sensitivity="asm10", maxiter=100)
    function verify(graph, msg="")
        if reference !== nothing
            for (name,path) ∈ graph.sequence
                seq = sequence(path)
                ref = reference[name]
                if seq != ref
                    badloci = Int[]
                    for i ∈ 1:min(length(seq),length(ref))
                        if seq[i] != ref[i]
                            push!(badloci, i)
                        end
                    end
                    left, right = max(badloci[1]-10, 1), min(badloci[1]+10, length(seq))
                    cumulative_lengths = cumsum([length(n.block, n) for n in path.node])

                    println("--> length:           ref($(length(ref))) <=> seq($(length(seq)))")
                    println("--> number of nodes:  $(length(path.node))")
                    println("--> path offset:      $(path.offset)")
                    println("--> |badloci|:        $(length(badloci))")
                    println("--> window:           $(left):$(badloci[1]):$(right)")
                    println("--> ref:              $(ref[left:right])") 
                    println("--> seq:              $(seq[left:right])") 

                    i = 1
                    while i < length(path.node)
                        if cumulative_lengths[i] > badloci[1]
                            break
                        else
                            i += 1
                        end
                    end

                    println("--> i(bad):           $(i)")
                    println("--> mutate:           $(path.node[i].block.mutate[path.node[i]])")
                    println("--> insert:           $(path.node[i].block.insert[path.node[i]])")
                    println("--> delete:           $(path.node[i].block.delete[path.node[i]])")

                    error("--> isolate '$name' incorrectly reconstructed")
                end
            end
        end

        graph
    end

    log("--> ordering")
    tree = ordering(compare, Gs...) |> balance
    log("--> tree: ", tree)

    meter = Progress(length(tree); desc="alignment progress", output=stderr)
    tips  = Dict{String,Graph}(collect(keys(G.sequence))[1] => G for G in Gs)

    log("--> aligning pairs")

    events = Channel{Bool}(10);
    result = Channel{Graph}(1);

    @spawn let
        while isopen(events)
            take!(events)
            next!(meter)
        end
    end

    G = nothing
    @sync for clade ∈ postorder(tree)
        @spawn begin
            if isleaf(clade)
                close(clade.graph)
                put!(clade.parent.graph, tips[clade.name])
                put!(events, true)
            else
                Gₗ = take!(clade.graph)
                Gᵣ = take!(clade.graph)
                close(clade.graph)

                G₀ = align_pair(Gₗ, Gᵣ, energy, minblock, verify, false; sensitivity=sensitivity)
                G₀ = align_self(G₀, energy, minblock, verify, false; sensitivity=sensitivity, maxiter=maxiter)

                put!(events, true)
                if clade.parent !== nothing
                    put!(clade.parent.graph, G₀)
                else
                    G = G₀
                end
            end
        end
    end
    close(events)

    return G
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
