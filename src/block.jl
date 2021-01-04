module Blocks

using Rematch, FStrings

using ..Utility: random_id, uncigar, wcpair, partition, Alignment
using ..Nodes

import ..Graphs: pair, reverse_complement

export SNPMap, InsMap, DelMap   # aux types
export Block 
export sequence, combine, swap! # operators

Maybe{T} = Union{T,Nothing}

# ------------------------------------------------------------------------
# Block data structure

# aliases
SNPMap = Dict{Int,UInt8}
InsMap = Dict{Tuple{Int,Int},Array{UInt8}} 
DelMap = Dict{Int,Int} 

mutable struct Block
    uuid::String
    sequence::Array{UInt8}
    gaps::Dict{Int,Int}
    mutate::Dict{Node{Block},SNPMap}
    insert::Dict{Node{Block},InsMap}
    delete::Dict{Node{Block},DelMap}
end

# ---------------------------
# constructors

# simple helpers
Block(sequence,gaps,mutate,insert,delete) = Block(random_id(),sequence,gaps,mutate,insert,delete)
Block(sequence) = Block(sequence,Dict{Int,Int}(),Dict{Node{Block},SNPMap}(),Dict{Node{Block},InsMap}(),Dict{Node{Block},DelMap}())
Block()         = Block(UInt8[])

translate(dict, δ) = Dict(key=>Dict(x+δ => v for (x,v) in val) for (key,val) in dict)
function translate!(dict, δ)
    for (key, val) in dict
        dict[key] = Dict(x+δ => v for (x,v) in val)
    end
end

# TODO: rename to concatenate?
# serial concatenate list of blocks
function Block(bs::Block...)
    @assert all([isolates(bs[1]) == isolates(b) for b in bs[2:end]])

    sequence = join([b.sequence for b in bs])

    gaps   = bs[1].gaps
    mutate = bs[1].mutate
    insert = bs[1].insert
    delete = bs[1].delete

    δ = length(bs[1])
    for b in bs[2:end]
        merge!(gaps,   translate(b.gaps,   δ))
        merge!(mutate, translate(b.mutate, δ))
        merge!(insert, translate(b.insert, δ))
        merge!(delete, translate(b.delete, δ))

        δ += length(b)
    end

    return Block(sequence,gaps,mutate,insert,delete)
end

# TODO: rename to slice?
# returns a subslice of block b
function Block(b::Block, slice)
    @assert slice.start >= 1 && slice.stop <= length(b)
    sequence = b.sequence[slice]

    select(dict,i) = translate(
                        Dict(node=>filter(p -> (first(p) >= i.start) && (first(p) <= i.stop), val) for (node,val) in dict), 
                     -i.start)

    gaps   = select(b.gaps,   slice)
    mutate = select(b.mutate, slice)
    insert = select(b.insert, slice)
    delete = select(b.delete, slice)

    return Block(sequence,gaps,mutate,insert,delete)
end

# ---------------------------
# operations

# simple operations
depth(b::Block) = length(b.mutate)
pair(b::Block)  = b.uuid => b

Base.show(io::IO, b::Block) = Base.show(io, (id=b.uuid, depth=depth(b)))

Base.length(b::Block)          = Base.length(b.sequence)
Base.length(b::Block, n::Node) = (length(b)
                               +((length(b.insert[n]) == 0) ? 0 : sum(length(i) for i in values(b.insert[n])))
                               -((length(b.delete[n]) == 0) ? 0 : sum(values(b.delete[n]))))

Locus = Union{
    NamedTuple{(:pos, :kind), Tuple{Int, Symbol}},
    NamedTuple{(:pos, :kind), Tuple{Tuple{Int,Int}, Symbol}},
}

islesser(a::Int, b::Int)                       = isless(a, b)
islesser(a::Tuple{Int,Int}, b::Int)            = isless(first(a), b)
islesser(a::Int, b::Tuple{Int,Int})            = isless(a, first(b))
islesser(a::Tuple{Int,Int}, b::Tuple{Int,Int}) = isless(a, b)

islesser(a::Locus, b::Locus) = islesser(a.pos, b.pos)

function allele_positions(b::Block, n::Node)
    keys(dict, sym) = [(pos=key, kind=sym) for key in Base.keys(dict)]
    return [keys(b.mutate[n],:snp); keys(b.insert[n],:ins); keys(b.delete[n],:del)]
end

# complex operations
function reverse_complement(b::Block)
    seq = reverse_complement(b.sequence)
    len = length(seq)

    revcmpl(dict::SNPMap) = Dict(len-locus+1:wcpair[nuc]  for (locus,nuc) in dict)
    revcmpl(dict::DelMap) = Dict(len-locus+1:del for (locus,del) in dict)
    revcmpl(dict::InsMap) = Dict((len-locus+1,b.gaps[locus]-off+1):reverse_complement(ins) for ((locus,off),ins) in dict)

    mutate = Dict(node => revcmpl(snp) for (node, snp) in b.mutate)
    insert = Dict(node => revcmpl(ins) for (node, ins) in b.insert)
    delete = Dict(node => revcmpl(del) for (node, del) in b.delete)
    gaps   = Dict(node => revcmpl(gap) for (node, gap) in b.gaps)

    return Block(seq,gaps,mutate,insert,delete)
end

function sequence(b::Block; gaps=false)
    if !gaps
        return copy(b.sequence)
    end
    
    len = length(b) + sum(values(b.gaps))
    seq = Array{UInt8}(undef, len)

    l, iₛ = 1, 1
    for r in sort(collect(keys(b.gaps)))
        len = r - l 
        seq[iₛ:iₛ+len] = b.seq[l:r]

        iₛ += len + 1
        len = b.gaps[r]
        seq[iₛ:iₛ+len] .= UInt8('-')

        l   = r + 1
        iₛ += len + 1
    end

    seq[iₛ:end] = b.sequence[l:end]

    return seq
end

# returns the sequence WITH mutations and indels applied to the consensus for a given tag 
function sequence(b::Block, node::Node{Block}; gaps=false)
    ref = sequence(b; gaps=gaps)
    len = gaps ? length(ref) : length(b, node)
    seq = Array{UInt8}('-'^len)

    loci = allele_positions(b, node)
    sort!(loci, lt=islesser)

    iᵣ, iₛ = 1, 1
    for l in loci
        δ = l.pos - iᵣ
        seq[iₛ:iₛ+δ-1] = ref[iᵣ:l.pos-1]
        iₛ += δ

        @match l.kind begin
            :snp => begin
                seq[iₛ] = b.mutate[node][l.pos]
                iₛ += 1
                iᵣ += δ + 1
            end
            :ins => begin
                # NOTE: insertions are indexed by the position they follow.
                #       since we stop 1 short, we finish here before continuing insertion.
                seq[iₛ] = ref[iᵣ]
                iₛ += 1

                ins = b.insert[node][l.pos]
                len = length(ins)

                seq[iₛ:iₛ+len] = ins

                iₛ += len + 1
                iᵣ += δ + 1
            end
            :del => begin
                # NOTE: deletions index the first position of the deletion. 
                #       this is the reason we stop 1 short above
                len = b.delete[node][l.pos]

                iₛ += gaps*(len+1)
                iᵣ  = l.pos + len + 1
            end
              _  => error("unrecognized locus kind")
        end
    end

    seq[iₛ:end] = ref[iᵣ:end]

    return seq
end

function Base.append!(b::Block, node::Node{Block}, snp::Maybe{SNPMap}, ins::Maybe{InsMap}, del::Maybe{DelMap})
    @assert node ∉ keys(b.mutate)
    @assert node ∉ keys(b.insert)
    @assert node ∉ keys(b.delete)

    if isnothing(snp)
        snp = SNPMap()
    end

    if isnothing(ins)
        ins = InsMap()
    end

    if isnothing(del)
        del = DelMap()
    end

    b.mutate[node] = snp
    b.insert[node] = ins
    b.delete[node] = del
end

function swap!(b::Block, oldkey::Node{Block}, newkey::Node{Block})
    b.mutate[newkey] = pop!(b.mutate, oldkey)
    b.insert[newkey] = pop!(b.insert, oldkey)
    b.delete[newkey] = pop!(b.delete, oldkey)
end

function swap!(b::Block, oldkey::Array{Node{Block}}, newkey::Node{Block})
    mutate = pop!(b.mutate, oldkey[1])
    insert = pop!(b.insert, oldkey[1])
    delete = pop!(b.delete, oldkey[1])

    for key in oldkey[2:end]
        merge!(mutate, pop!(b.mutate, key))
        merge!(insert, pop!(b.insert, key))
        merge!(delete, pop!(b.delete, key))
    end

    b.mutate[newkey] = mutate
    b.insert[newkey] = insert
    b.delete[newkey] = delete 
end

function reconsensus!(b::Block)
    #nop
end

function combine(qry::Block, ref::Block, aln::Alignment; maxgap=500)
    sequences,intervals,mutations,inserts,deletes = partition(
                                                         uncigar(aln.cigar),
                                                         qry.sequence,
                                                         ref.sequence,
                                                         maxgap=maxgap
                                                    )

    blocks = NamedTuple{(:block,:kind),Tuple{Block,Symbol}}[]

    for (seq,pos,snp,ins,del) in zip(sequences,intervals,mutations,inserts,deletes)
        @match (pos.qry, pos.ref) begin
            ( nothing, rₓ )  => push!(blocks, (block=Block(ref, rₓ), kind=:ref))
            ( qₓ , nothing ) => push!(blocks, (block=Block(qry, qₓ), kind=:qry))
            ( qₓ , rₓ )      => begin
                @assert !isnothing(snp)
                @assert !isnothing(ins)
                @assert !isnothing(del)

                # slice both blocks
                r = Block(ref, rₓ)
                q = Block(qry, qₓ)

                # apply global snp and indels to all query sequences
                # XXX: do we have to worry about overlapping insertions/deletions?
                for node in keys(q.mutate)
                    merge!(q.mutate[node],snp)
                    merge!(q.insert[node],ins)
                    merge!(q.delete[node],del)
                end

                gap = Dict(first(key)=>length(val) for (key,val) in ins)
                new = Block(seq,gap,snp,ins,del)
                reconsensus!(new)

                push!(blocks, (block=new, kind=:all))
            end
        end
    end

    return blocks
end

# ------------------------------------------------------------------------
# main point of entry for testing

using Random, Distributions, StatsBase

function generate_alignment(;len=100,num=10,μ=(snp=1e-2,ins=0.0,del=1e-2),Δ=5)
    ref = Array{UInt8}(random_id(;len=len, alphabet=['A','C','G','T']))
    aln = zeros(UInt8, num, len)

    map = (
        snp = Array{SNPMap}(undef,num),
        ins = Array{InsMap}(undef,num),
        del = Array{DelMap}(undef,num),
    )
    ρ  = (
        snp = Poisson(μ.snp*len),
        ins = Poisson(μ.ins*len),
        del = Poisson(μ.del*len),
    )
    n = (
        snp = rand(ρ.snp, num),
        ins = rand(ρ.ins, num),
        del = rand(ρ.del, num),
    )

    for i in 1:num
        aln[i,:] = ref

        index = collect(1:len)

        # deletions
        loci = sample(index, n.del[i]; replace=false)
        dels = min.(len .- loci, sample(1:Δ, n.del[i]))
        for (locus, del) in zip(loci,dels)
            aln[i,locus:locus+del] .= UInt8('-')
        end

        map.del[i] = DelMap(zip(loci,dels))
        
        # single nucleotide polymorphisms
        loci = sample(1:len, n.snp[i]; replace=false)
        snps = sample(UInt8['A','C','G','T'], n.snp[i])

        redo = findall(snps .== ref[loci])
        while length(redo) >= 1
            snps[redo] = sample(UInt8['A','C','G','T'], length(redo))
            redo = findall(snps .== ref[loci])
        end

        for (locus,snp) in zip(loci,snps)
            aln[i,locus] = snp
        end

        map.snp[i] = SNPMap(zip(loci,snps))
    end

    return ref, aln, map
end

function test()
    Random.seed!(0)
    ref, aln, map = generate_alignment()
    to_char(d::Dict{Int,UInt8}) = Dict{Int,Char}(k=>Char(v) for (k,v) in d)

    blk  = Block(ref)
    node = [Node{Block}(blk,true) for i in 1:size(aln,1)]
    for i in 1:size(aln,1)
        append!(blk, node[i], map.snp[i], nothing, map.del[i])
    end

    ok  = true
    pos = join([f"{i:02d}" for i in 1:10:100], ' '^8)
    for i in 1:size(aln,1)
        seq  = sequence(blk,node[i];gaps=true)
        good = aln[i,:] .== seq
        if !all(good)
            ok = false

            err        = copy(seq)
            err[good] .= ' '

            println(f"failure on row {i}")
            println("Loci: ", pos)
            println("True: ", String(copy(aln[i,:])))
            println("Estd: ", String(copy(seq)))
            println("Diff: ", String(err))
            println("SNPs: ", to_char(map.snp[i]))
            println("Dels: ", map.del[i])
        end
    end

    if ok
        println("worked!")
    else
        println("failed: check diagnostics!")
    end
end

end
