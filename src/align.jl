module Align

using Rematch, Dates
using LinearAlgebra
using ProgressMeter

import Base.Threads.@spawn

using ..Utility: read_paf, enforce_cutoff!, reverse_complement
using ..Blocks
using ..Paths: replace!
using ..Nodes
using ..Graphs

# Pool belongs only to the alignment module!
include("pool.jl")
using .Pool

export align

# ------------------------------------------------------------------------
# global variables
# TODO: move to a better location

fifos       = pool(2) #pool(3*Threads.nthreads()) NOTE: uncomment when you want concurrency
getio()     = take!(fifos)
hasio()     = isready(fifos)
putio(fifo) = put!(fifos, fifo)
finalize()  = shutdown(fifos)
atexit(finalize)

# TODO: generalize to n > 2
# NOTE: this is to ensure no deadlock for singleton pulls
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
    println(stderr, now(), " ", join(msg, " ")...)
    flush(stderr)
end

# command line execution
function execute(cmd::Cmd; now=true)
    out = Pipe()
    err = Pipe()

    proc = run(pipeline(cmd, stdout=out, stderr=err); wait=now)

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
    return execute(`minimap2 -x asm10 -m 10 -n 1 -s 30 -D -c $ref $qry`; now=false)
end

# ------------------------------------------------------------------------
# guide tree for order of pairwise comparison for multiple alignments

# TODO: distance?
mutable struct Clade
    name   :: String
    parent :: Union{Clade,Nothing}
    left   :: Union{Clade,Nothing}
    right  :: Union{Clade,Nothing}
    graph  :: Channel{Graph}
end

# ---------------------------
# constructors

Clade()     = Clade("",nothing,nothing,nothing,Channel{Graph}(1))
Clade(name) = Clade(name,nothing,nothing,nothing,Channel{Graph}(1))
Clade(left::Clade, right::Clade) = Clade("",nothing,left,right,Channel{Graph}(1))

function Clade(distance, names; algo=:nj)
    @match algo begin
        :nj => return nj(distance, names)
        _ => error("unrecognized algorithm $(algo)")
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

function Base.length(node::Clade)
    if isleaf(node)
        return 1
    else
        return 1 + Base.length(node.left) + Base.length(node.right)
    end
end

# ------------------------------------------------------------------------
# ordering functions

# TODO: assumes the input graphs are singletons! generalize
function ordering(Gs...; compare=mash)
    fifo = getio()

    task = @async compare(path(fifo))

    @async begin
        open(fifo,"w") do io
            for G ∈ Gs
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

function write(fifo, G::Graph)
    io = open(fifo, "w")

    @label write
    sleep(1e-3) # NOTE: hack to allow for minimap2 to open file
    try
        marshal(io, G)
    catch e
        log("ERROR: $(e)")
        @goto write
        # error(e)
    finally
        close(io)
    end
end

function do_align(G₁::Graph, G₂::Graph, io₁, io₂, energy::Function)
    # NOTE: minimap2 opens up file descriptors in order!
    #       must process 2 before 1 otherwise we deadlock
    
    cmd = minimap2(path(io₁), path(io₂))
    @async begin
        write(io₂, G₂) # ref
        write(io₁, G₁) # qry
    end

    out = IOBuffer(fetch(cmd.out)) # NOTE: blocks until minimap finishes

    hits = collect(read_paf(out))
    sort!(hits; by=energy)

    close(out)

    return hits
end

function align_kernel(hits, energy, minblock, skip, blocks!, replace!)
    blocks = Dict{String,Block}()
    ok = false
    for hit in hits
        energy(hit) >= 0 && break
        skip(hit) && continue

        ok   = true
        qry₀, ref₀ = blocks!(hit)

        hit.qry.seq = qry₀.sequence
        hit.ref.seq = ref₀.sequence
        enforce_cutoff!(hit, minblock)

        # log(hit)

        blks, strand = combine(qry₀, ref₀, hit; minblock=minblock)

        qrys = map(b -> b.block, filter(b -> b.kind != :ref, blks))
        refs = map(b -> b.block, filter(b -> b.kind != :qry, blks))

        replace!((qry=qry₀, ref=ref₀), (qry=qrys, ref=refs), strand)

        for blk in map(b->b.block, blks)
            blocks[blk.uuid] = blk
        end
    end

    return blocks, ok
end

function align_self(G₁::Graph, io₁, io₂, energy::Function, minblock::Int, verify::Function; maxiter=100)
    G₀ = G₁
    ok = true

    niter = 0
    while ok
        ok   = false
        hits = do_align(G₀, G₀, io₁, io₂, energy)
        
        skip  = (hit) -> (hit.qry.name == hit.ref.name) || (!(hit.qry.name in keys(G₀.block)) || !(hit.ref.name in keys(G₀.block)))
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

        blocks, ok = align_kernel(hits, energy, minblock, skip, block, replace)
        
        merge!(blocks, G₀.block)

        if ok && niter < maxiter
            G₀ = Graph(
                blocks,
                G₀.sequence,
            )
            detransitive!(G₀)
            niter += 1
        end
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

function align_pair(G₁::Graph, G₂::Graph, io₁, io₂, energy::Function, minblock::Int)
    hits = do_align(G₁, G₂, io₁, io₂, energy)

    skip  = (hit) -> !(hit.qry.name in keys(G₁.block)) || !(hit.ref.name in keys(G₂.block))
    block = (hit) -> (
        qry = pop!(G₁.block, hit.qry.name),
        ref = pop!(G₂.block, hit.ref.name),
    )
    replace = (old, new, orientation) -> let
        # START DEBUG
        compare(old.qry, new.qry, orientation)
        compare(old.ref, new.ref, true)
        # END DEBUG

        for path in values(G₁.sequence)
            replace!(path, old.qry, new.qry, orientation)
        end

        for path in values(G₂.sequence)
            replace!(path, old.ref, new.ref, true)
        end
    end

    blocks, _ = align_kernel(hits, energy, minblock, skip, block, replace)
    sequence  = merge(G₁.sequence, G₂.sequence)

    # XXX: worry about uuid collision?
    merge!(blocks, G₁.block)
    merge!(blocks, G₂.block)

    G = Graph(
        blocks, 
        sequence,
    )
    detransitive!(G)

    return G
end

# TODO: the associative array is a bit hacky...
#       can we push it directly into the channel?
function align(Gs::Graph...; energy=(hit)->(-Inf), minblock=100, reference=nothing)
    function verify(graph, msg="")
        if reference !== nothing
            log(msg)

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
    tree = ordering(Gs...)
    log("--> tree: ", tree)

    meter = Progress(length(tree); desc="alignment progress", output=stderr)
    function kernel(clade)
        Gₗ = take!(clade.left.graph)
        Gᵣ = take!(clade.right.graph)

        io₁, io₂ = getios()

        G₀ = align_pair(Gₗ, Gᵣ, io₁, io₂, energy, minblock)
        G₀ = align_self(G₀, io₁, io₂, energy, minblock, verify)

        putio(io₁)
        putio(io₂)

        next!(meter)
        put!(clade.graph, G₀)
    end

    # sequences on tips of tree
    tips = Dict{String,Graph}(collect(keys(G.sequence))[1] => G for G in Gs)

    log("--> aligning pairs")
    for clade ∈ postorder(tree)
        if isleaf(clade)
            put!(clade.graph, tips[clade.name])
            next!(meter)
            close(clade.graph)
        else
            kernel(clade)
            close(clade.graph)
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
