module Blocks

using Rematch

import Base:
    show, length, append!, keys

# internal modules
using ..Intervals
using ..Nodes
using ..Utility: 
    random_id, contiguous_trues,
    uncigar, wcpair, Alignment

import ..Graphs:
    pair, reverse_complement, 
    sequence, sequence!

# exports
export SNPMap, InsMap, DelMap # aux types
export Block 
export combine, swap!, check  # operators

# ------------------------------------------------------------------------
# utility types

Maybe{T} = Union{T,Nothing}

# aliases
const SNPMap = Dict{Int,UInt8}
const InsMap = Dict{Tuple{Int,Int},Array{UInt8}} 
const DelMap = Dict{Int,Int} 

show(io::IO, m::SNPMap) = show(io, [ k => Char(v) for (k,v) in m ])
show(io::IO, m::InsMap) = show(io, [ k => String(Base.copy(v)) for (k,v) in m ])

# ------------------------------------------------------------------------
# utility functions

mutable struct Pos
    start::Int
    stop::Int
end

Base.to_index(x::Pos) = x.start:x.stop
advance!(x::Pos)      = x.start=x.stop
copy(x::Pos)          = Pos(x.start,x.stop)

mutable struct PairPos
    qry::Maybe{Pos}
    ref::Maybe{Pos}
end

# TODO: relax hardcoded reliance on cigar suffixes. make symbols instead
const PosPair    = NamedTuple{(:qry, :ref), Tuple{Maybe{Pos}, Maybe{Pos}}} 
const BlockAlign = NamedTuple{(:range, :segment), Tuple{PosPair, PosPair[]}}

function partition(alignment; maxgap=500)
    qry, ref = alignment.qry.seq, alignment.ref.seq

    qryₓ = Pos(1,1)
    refₓ = Pos(1,1)

    block   = BlockAlign[]
    segment = PosPair[]  # segments of current block being constructed

    # ----------------------------
    # internal operators
    
    function finalize_block!()
        length(segment) == 0 && @goto advance

        push!(block, (
            range   = (
                qry = Pos(qryₓ.start,qryₓ.stop-1), 
                ref = Pos(refₓ.start,refₓ.stop-1)
            ),
            segment = segment
        ))

        segment = PosPair[]
        
        @label advance
        advance!(qryₓ)
        advance!(refₓ)
    end

    function qry_block!(pos)
        push!(block, (
            range   = (
                qry = pos,
                ref = nothing
            ),
            segment = PosPair[]
         ))
    end

    function ref_block!(pos)
        push!(block, (
            range   = (
                qry = nothing,
                ref = pos,
            ),
            segment = PosPair[]
         ))
    end

    # ----------------------------
    # see if blocks have a leading unmatched block

    if alignment.qry.start > 1
        qry_block!(Pos(1,alignment.qry.start-1))
        qryₓ = Pos(alignment.qry.start,alignment.qry.start)
    end

    if alignment.ref.start > 1
        ref_block!(Pos(1, alignment.ref.start-1))
        refₓ = Pos(alignment.ref.start,alignment.ref.start)
    end
    
    # ----------------------------
    # parse cigar within region of overlap
    
    for (len, type) ∈ uncigar(alignment.cigar)
        @match type begin
        'S' || 'H' => begin
            # XXX:  treat soft clips differently?
            # TODO: implement
            error("need to implement soft/hard clipping")
        end
        'M' => begin
            x = Pos(refₓ.stop, refₓ.stop+len-1)
            y = Pos(qryₓ.stop, qryₓ.stop+len-1)

            push!(segment, (ref=x, qry=y))

            qryₓ.stop += len
            refₓ.stop += len
        end
        'D' => begin
            if len >= maxgap
                finalize_block!()

                ref_block!(Pos(refₓ.start,refₓ.stop+len-1))

                refₓ.stop += len
                advance!(refₓ)
            else
                push!(segment, (ref=Pos(refₓ.stop, refₓ.stop+len-1), qry=nothing))
                refₓ.stop += len
            end
        end
        'I' => begin
            if len >= maxgap
                finalize_block!()

                qry_block!(Pos(qryₓ.start,qryₓ.stop+len-1))

                qryₓ.stop += len
                advance!(qryₓ)
            else
                push!(segment, (ref=nothing, qry=Pos(qryₓ.stop,qryₓ.stop+len-1)))
                qryₓ.stop += len
            end
        end
         _  => error("unrecognized cigar string suffix")
        end
    end

    finalize_block!()

    # ----------------------------
    # see if blocks have a trailing unmatched block

    if alignment.qry.stop < alignment.qry.length
        qry_block!(Pos(alignment.qry.stop,alignment.qry.length))
    end

    if alignment.ref.stop < alignment.ref.length
        ref_block!(Pos(alignment.ref.stop,alignment.ref.length))
    end

    return block
end

# ------------------------------------------------------------------------
# Block data structure

mutable struct Block
    uuid     :: String
    sequence :: Array{UInt8}
    gaps     :: Dict{Int,Int}
    mutate   :: Dict{Node{Block},SNPMap}
    insert   :: Dict{Node{Block},InsMap}
    delete   :: Dict{Node{Block},DelMap}
end

function show(io::IO, m::Dict{Node{Block}, T}) where T <: Union{SNPMap, InsMap, DelMap}
    print(io, "{\n")
    for (k,v) in m
        print(io, "\t", k, " => {")
        show(io, v)
        print(io, "}\n")
    end
    print(io, "}\n")
end

# ---------------------------
# constructors

# simple helpers
Block(sequence,gaps,mutate,insert,delete) = Block(random_id(),sequence,gaps,mutate,insert,delete)
Block(sequence) = Block(sequence,Dict{Int,Int}(),Dict{Node{Block},SNPMap}(),Dict{Node{Block},InsMap}(),Dict{Node{Block},DelMap}())
Block()         = Block(UInt8[])

translate(d::Dict{Int,Int}, δ) = Dict(x+δ => v for (x,v) ∈ d) # gaps
translate(d::Dict{Node{Block},InsMap}, δ) = Dict(n => Dict((x+δ,Δ) => v for ((x,Δ),v) ∈ val) for (n,val) ∈ d) # insertions 
translate(dict, δ) = Dict(key=>Dict(x+δ => v for (x,v) in val) for (key,val) in dict)

# TODO: rename to concatenate?
# serial concatenate list of blocks
function Block(bs::Block...)
    sequence = vcat((b.sequence for b in bs)...)

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
    if (slice.start == 1 && slice.stop == length(b))
        Block(b.sequence,b.gaps,b.mutate,b.insert,b.delete)
    end
    @assert slice.start >= 1 && slice.stop <= length(b)

    sequence = b.sequence[slice]

    select(dict,i) = translate(
        Dict(node => filter(
            (p) -> (i.start ≤ first(first(p)) ≤ i.stop), 
            subdict) for (node, subdict) ∈ dict
        ), 
        -i.start+1
    )

    gaps = Dict(x-slice.start+1 => δ for (x,δ) ∈ b.gaps if slice.start ≤ x ≤ slice.stop)

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

show(io::IO, b::Block) = show(io, (id=b.uuid, depth=depth(b)))

length(b::Block) = length(b.sequence)
length(b::Block, n::Node) = (length(b)
                          +((length(b.insert[n]) == 0) ? 0 : sum(length(i) for i in values(b.insert[n])))
                          -((length(b.delete[n]) == 0) ? 0 : sum(values(b.delete[n]))))

keys(b::Block) = keys(b.mutate)

# internal structure to allow us to sort all allelic types
Locus = Union{
    NamedTuple{(:pos, :kind), Tuple{Int, Symbol}},
    NamedTuple{(:pos, :kind), Tuple{Tuple{Int,Int}, Symbol}},
}

islesser(a::Int, b::Int)                       = isless(a, b)
islesser(a::Tuple{Int,Int}, b::Int)            = isless(first(a), b)
islesser(a::Int, b::Tuple{Int,Int})            = isless(a, first(b)) || a == first(b)
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
    !gaps && return b.sequence
    
    len = length(b) + sum(values(b.gaps))
    seq = Array{UInt8}(undef, len)

    l, iₛ = 1, 1
    for r in sort(collect(keys(b.gaps)))
        len = r - l
        seq[iₛ:iₛ+len] = b.sequence[l:r]

        iₛ += len + 1
        len = b.gaps[r]
        seq[iₛ:iₛ+len-1] .= UInt8('-')

        l   = r + 1
        iₛ += len
    end

    seq[iₛ:end] = b.sequence[l:end]

    return seq
end

function sequence_gaps!(seq, b::Block, node::Node{Block})
    ref = sequence(b; gaps=true)
    @assert length(seq) == length(ref)

    loci = allele_positions(b, node) 
    sort!(loci, lt=islesser)

    Ξ(x) = x + reduce(+,(δ for (l,δ) in b.gaps if l < x); init=0)

    @show b.uuid
    @show b.gaps
    @show length(b.sequence), length(seq), length(ref), length(b, node)
    @show b.mutate[node]
    @show b.insert[node]
    @show b.delete[node]

    for l in loci
        @match l.kind begin
            :snp => begin
                x         = l.pos
                seq[Ξ(x)] = b.mutate[node][x]
            end
            :ins => begin
                ins = b.insert[node][l.pos]
                len = length(ins)

                x = Ξ(l.pos[1]) # NOTE: insertion occurs 1 nt AFTER the key position
                δ = l.pos[2]

                seq[x+δ+1:x+len+δ] = ins
            end
            :del => begin
                len = b.delete[node][l.pos]
                x   = Ξ(l.pos )

                @show l, x, len, length(seq)

                seq[x:x+len-1] .= UInt8('-')
            end
              _  => error("unrecognized locus kind")
        end
    end

    return seq
end

function sequence_gaps(b::Block, node::Node{Block})
    len = length(b) + sum(values(b.gaps)) # TODO: make alignment_length function?
    seq = Array{UInt8}(undef, len)

    sequence_gaps!(seq, b, node)

    return seq
end

# returns the sequence WITH mutations and indels applied to the consensus for a given tag 
function sequence!(seq, b::Block, node::Node{Block}; gaps=false)
    gaps && return sequence_gaps!(seq, b, node)

    @assert length(seq) == length(b, node)

    ref = sequence(b; gaps=false)

    pos  = (l) -> isa(l.pos, Tuple) ? l.pos[1] : l.pos # dispatch over different key types
    loci = allele_positions(b, node)
    sort!(loci, lt=islesser)

    iᵣ, iₛ = 1, 1
    for l in loci
        if (δ = pos(l) - iᵣ) >= 0
            seq[iₛ:iₛ+δ-1] = ref[iᵣ:pos(l)-1]
            iₛ += δ
        end

        @match l.kind begin
            :snp => begin
                seq[iₛ] = b.mutate[node][l.pos]
                iₛ += 1
                iᵣ += δ + 1
            end
            :ins => begin
                # NOTE: insertions are indexed by the position they follow.
                #       since we stop 1 short, we finish here before continuing insertion.
                if δ >= 0
                    seq[iₛ] = ref[pos(l)]
                    iₛ += 1
                end

                ins = b.insert[node][l.pos]
                len = length(ins)

                seq[iₛ:iₛ+len-1] = ins

                iₛ += len
                iᵣ  = pos(l) + 1
            end
            :del => begin
                # NOTE: deletions index the first position of the deletion. 
                #       this is the reason we stop 1 short above
                iᵣ = l.pos + b.delete[node][l.pos]
            end
              _  => error("unrecognized locus kind")
        end
    end

    seq[iₛ:end] = ref[iᵣ:end]

    return seq
end

function sequence(b::Block, node::Node{Block}; gaps=false)
    seq = gaps ? sequence(b; gaps=true) : Array{UInt8}('-'^length(b, node))
    sequence!(seq, b, node; gaps=gaps)
    return seq
end

function append!(b::Block, node::Node{Block}, snp::Maybe{SNPMap}, ins::Maybe{InsMap}, del::Maybe{DelMap})
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
    # NOTE: no point to compute this for blocks with 1 or 2 individuals
    depth(b) <= 2 && return false 

    # NOTE: we can't assume that keys(b) will return the same order on subsequent calls
    #       thus we collect into array here for a static ordering of the nodes
    nodes = collect(keys(b))

    ref = sequence(b; gaps=true)
    aln = Array{UInt8}(undef, length(ref), depth(b))
    for (i,node) in enumerate(nodes)
        aln[:,i] = ref
        sequence!(view(aln,:,i), b, node; gaps=true)
    end

    consensus = [mode(view(aln,i,:)) for i in 1:size(aln,1)]
    if all(consensus .== ref) # hot path: if consensus sequence did not change, abort!
        return false
    end

    isdiff = (aln .!= consensus)
    refdel = (consensus .== UInt8('-'))
    alndel = (aln .== UInt8('-'))

    δ = (
        snp = isdiff .& .!refdel .& .!alndel,
        del = isdiff .& .!refdel .&   alndel,
        ins = isdiff .&   refdel .& .!alndel,
    )

    coord   = cumsum(.!refdel)

    refgaps = contiguous_trues(refdel)
    b.gaps  = Dict{Int, Int}(coord[gap.lo] => length(gap) for gap in refgaps)
    
    b.mutate = Dict{Node{Block},SNPMap}( 
            node => SNPMap(
                      coord[l] => aln[l,i] 
                for l in findall(δ.snp[:,i])
            )
        for (i,node) in enumerate(nodes)
    )

    b.delete = Dict{Node{Block},DelMap}( 
            node => DelMap(
                      coord[del.lo] => length(del)
                for del in contiguous_trues(δ.del[:,i])
             )
        for (i,node) in enumerate(nodes)
    )

    Δ(I) = (R = containing(refgaps, I)) == nothing ? 0 : I.lo - R.lo
    b.insert = Dict{Node{Block},InsMap}( 
            node => InsMap(
                      (coord[ins.lo],Δ(ins)) => aln[ins,i] 
                for ins in contiguous_trues(δ.ins[:,i])
             )
        for (i,node) in enumerate(nodes)
    )

    b.sequence = consensus[.!refdel]
    @show b.uuid
    @show length(b.sequence), length(consensus)
    @show b.mutate
    @show b.insert
    @show b.delete

    @assert all(all(k ≤ length(b.sequence) for k in keys(d)) for d in values(b.mutate)) 
    @assert all(all(k ≤ length(b.sequence) for k in keys(d)) for d in values(b.delete)) 
    @assert all(all(k[1] ≤ length(b.sequence) for k in keys(d)) for d in values(b.insert)) 

    return true
end

function combine(qry::Block, ref::Block, aln::Alignment; maxgap=500)
    # NOTE: this will enforce that indels are less than maxgap!
    # TODO: rename partition function
    sequences,intervals,mutations,inserts,deletes = partition(aln; maxgap=maxgap)

    blocks = NamedTuple{(:block,:kind),Tuple{Block,Symbol}}[]
    for (seq,pos,snp,ins,del) in zip(sequences,intervals,mutations,inserts,deletes)
        @match (pos.qry, pos.ref) begin
            ( nothing, rₓ )  => begin
                push!(blocks, (block=Block(ref, rₓ), kind=:ref))
            end
            ( qₓ , nothing ) => begin
                push!(blocks, (block=Block(qry, qₓ), kind=:qry))
            end
            ( qₓ , rₓ )      => begin
                @assert !isnothing(snp)
                @assert !isnothing(ins)
                @assert !isnothing(del)

                # slice both blocks to window of overlap
                @show ref, qry
                r = Block(ref, rₓ)
                q = Block(qry, qₓ)

                @assert all(r.sequence .== seq) # NOTE: we can most likely get rid of sequence

                @show snp
                @show ins
                @show del
                @show q.mutate
                @show q.insert
                @show q.delete

                rereference!(q, r, (mutate=snp, insert=ins, delete=del))

                gaps = Dict(first(key)=>length(val) for (key,val) in ins)
                new  = Block(
                    seq,
                    merge(r.gaps,q.gaps,gaps),
                    merge(r.mutate,q.mutate),
                    merge(r.insert,q.insert),
                    merge(r.delete,q.delete),
                )

                @show rₓ, qₓ
                @show new.uuid
                @show length(new.sequence), length(r.sequence), length(q.sequence)

                @assert all(all(k ≤ length(new.sequence) for k in keys(d)) for d in values(new.mutate)) 
                @assert all(all(k ≤ length(new.sequence) for k in keys(d)) for d in values(new.delete)) 
                @assert all(all(k[1] ≤ length(new.sequence) for k in keys(d)) for d in values(new.insert)) 

                reconsensus!(new)

                push!(blocks, (block=new, kind=:all))
            end
        end
    end

    return blocks
end

# XXX: we have to worry about overlapping insertions/deletions
#      consider the following case:
#
#      Ref: *****
#      Qry: *****A
#
#      we would store the last A as a global insertion for all qrys
#      however, what if the 'A' itself is deleted within a minority of contained sequences?
#      the same would be true if Ref had the extra A and a few qrys contained this insertion
#      thus we need some way for merged insertions and deletions to annihilate
#      furthermore, it does not need to occur at the end, can occur anywhere in sequence
#
#      this suggests a general algorithm: 
#      we need to iterate through all contained sequences in qry
#      -> ask if any portion of a global insertion within ins is deleted. 
#      ---> if so, make adjustments
#      -> check if any global deletion has an insertion within.
#      -> check if any global insertion is modified by a snp
#
# NOTE: allele dictionaries are referenced to position on new consensus sequence
function rereference!(qry::Block, ref::Block, node::Node, allele)
    #=
    # new -> old coordinate map
    Ξₒ(xₙ) = (xₙ
           + sum(length(v) for (x,v) in allele.insert if first(x) < xₙ)
           - sum(v for (x,v) in allele.delete if x < xₙ))
    # old -> new coordinate map
    Ξₙ(xₒ) = (xₒ
           - sum(length(v) for (xₙ,v) in allele.insert if Ξₒ(first(xₙ)) < xₒ)
           + sum(v for (xₙ,v) in allele.delete if Ξₒ(xₙ) < xₒ))
    =#

    # seq = sequence(b, node; gaps=true)

    merge!(qry.mutate[node], allele.mutate)
    merge!(qry.insert[node], allele.insert)
    merge!(qry.delete[node], allele.delete)
end

rereference!(qry::Block, ref::Block, allele) = for node in keys(b) 
    rereference!(qry, ref, node, allele)
end

function check(b::Block)
    @show b.uuid
    @show b.mutate
    @show b.insert
    @show b.delete
    @assert all( n.block == b for n ∈ keys(b) )
    @assert all( n.block == b for n ∈ keys(b) )
    @assert all( n.block == b for n ∈ keys(b) )
end

# ------------------------------------------------------------------------
# main point of entry for testing

using Random, Distributions, StatsBase

function generate_alignment(;len=100,num=10,μ=(snp=1e-2,ins=1e-2,del=1e-2),Δ=5)
    ref = Array{UInt8}(random_id(;len=len, alphabet=['A','C','G','T']))
    aln = zeros(UInt8, num, len)

    map = (
        snp = Array{SNPMap}(undef,num),
        ins = Array{InsMap}(undef,num),
        del = Array{DelMap}(undef,num),
    )
    ρ = (
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
    end

    # random insertions
    # NOTE: this is the inverse operation as a deletion.
    #       perform operation as a collective.
    inserts = Array{IntervalSet{Int}}(undef, num)

    # first collect all insertion intervals
    for i in 1:num
        inserts[i] = IntervalSet(1, len+1)

        for j in 1:n.ins[i]
            @label getinterval
            start = sample(1:len)
            delta = len-start+1
            stop  = start + min(delta, sample(1:Δ))

            insert = Interval(start, stop)

            if !isdisjoint(inserts[i], insert)
                @goto getinterval # XXX: potential infinite loop
            end

            inserts[i] = inserts[i] ∪ insert
        end
    end

    allinserts = reduce(∪, inserts)

    δ = 1 
    gaps = [begin 
        x  = (I.lo-δ, length(I)) 
        δ += length(I)
        x
    end for I in allinserts]

    for (i, insert) in enumerate(inserts)
        keys = Array{Tuple{Int,Int}}(undef, length(insert))
        vals = Array{Array{UInt8}}(undef, length(insert))
        for (n, a) in enumerate(insert)
            for (j, b) in enumerate(allinserts)
                if a ⊆ b
                    keys[n] = (gaps[j][1], a.lo - b.lo)
                    vals[n] = ref[a]
                    @goto outer
                end
            end
            error("failed to find containing interval!")
            @label outer
        end

        map.ins[i] = InsMap(zip(keys,vals))

        # delete non-overlapping regions
        for j in allinserts \ insert
            aln[i,j] .= UInt8('-')
        end
    end

    idx = collect(1:len)[~allinserts]
    ref = ref[~allinserts]

    for i in 1:num
        index = collect(1:length(idx))
        deleteat!(index, findall(aln[i,idx] .== UInt8('-')))

        # random deletions
        # NOTE: care must be taken to ensure that they don't overlap or merge
        loci = Array{Int}(undef, n.del[i])
        dels = Array{Int}(undef, n.del[i])

        for j in 1:n.del[i]
            @label tryagain
            loci[j] = sample(index)

            while aln[i,max(1, idx[loci[j]]-1)] == UInt8('-')
                loci[j] = sample(index)
            end

            x = idx[loci[j]]

            offset = findfirst(aln[i,x:end] .== UInt8('-'))
            maxgap = isnothing(offset) ? (len-x+1) : (offset-1)

            dels[j] = min(maxgap, sample(1:Δ))

            # XXX: this is a hack to ensure deletions and insertions don't overlap
            if !all(item ∈ idx for item in x:x+dels[j]-1)
                @goto tryagain
            end

            aln[i,x:(x+dels[j]-1)] .= UInt8('-')
            filter!(i->i ∉ loci[j]:(loci[j]+dels[j]-1), index)
        end

        map.del[i] = DelMap(zip(loci,dels))
        
        # random single nucleotide polymorphisms
        # NOTE: we exclude the deleted regions
        loci = sample(index, n.snp[i]; replace=false)
        snps = sample(UInt8['A','C','G','T'], n.snp[i])
        redo = findall(ref[loci] .== snps)

        while length(redo) >= 1
            snps[redo] = sample(UInt8['A','C','G','T'], length(redo))
            redo = findall(ref[loci] .== snps)
        end

        for (locus,snp) in zip(loci,snps)
            aln[i,idx[locus]] = snp
        end

        map.snp[i] = SNPMap(zip(loci,snps))
    end

    return ref, aln, Dict(gaps), map
end

function verify(blk, node, aln)
    local pos = join(["$(i)" for i in 1:10:101], ' '^8)
    local tic = join(["|" for i in 1:10:101], '.'^9)

    ok = true
    for i in 1:size(aln,1)
        seq  = sequence(blk,node[i];gaps=true)
        if size(aln,2) != length(seq)
            println("failure on row $(i), node $(node[i])")
            println("incorrect size!")
            ok = false
            break
        end

        good = aln[i,:] .== seq
        if !all(good)
            ok = false

            err        = copy(seq)
            err[good] .= ' '

            println("failure on row $(i), node $(node[i])")
            println("Loci: ", pos)
            println("      ", tic)
            println("Ref:  ", String(copy(sequence(blk; gaps=true))))
            println("True: ", String(copy(aln[i,:])))
            println("Estd: ", String(copy(seq)))
            println("Diff: ", String(err))
            println("SNPs: ", blk.mutate[node[i]])
            println("Dels: ", blk.delete[node[i]])
            println("Ints: ", blk.insert[node[i]])
            break
        end
        seq  = sequence(blk,node[i];gaps=false)
    end

    return ok
end

function test()
    ref, aln, gap, map = generate_alignment()

    blk = Block(ref)
    blk.gaps = gap

    node = [Node{Block}(blk,true) for i in 1:size(aln,1)]
    for i in 1:size(aln,1)
        append!(blk, node[i], map.snp[i], map.ins[i], map.del[i])
    end

    ok = verify(blk, node, aln)
    if !ok
        error("failure to initialize block correctly")
    end

    reconsensus!(blk)

    ok = verify(blk, node, aln)
    if !ok
        error("failure to reconsensus block correctly")
    end

    return ok 
end

end
