module Utility

using FStrings, Match
using StatsBase

# NOTE: for debugging/benchmarking
using BenchmarkTools, Infiltrator
using Profile, ProfileView

import Base.Threads.@spawn

import ..Graphs: reverse_complement

export random_id, log
export Alignment
export enforce_cutoff, cigar, uncigar, homologous

export read_fasta, name
export read_paf

Maybe{T} = Union{Nothing,T}

# ------------------------------------------------------------------------
# random functions

# random string of fixed length
function random_id(;len=10, alphabet=UInt8[])
    if length(alphabet) == 0
        alphabet = ['A','B','C','D','E','F','G','H','I','J','K','L','M',
                    'N','O','P','Q','R','S','T','U','V','W','X','Y','Z']
    end
    return join(sample(alphabet, len))
end

# ------------------------------------------------------------------------
# banded alignment data structure

struct Score <: AbstractArray{Float64,2}
    rows::Int
    cols::Int
    band::NamedTuple{(:lower, :upper)}
    data::Array{Float64}
    offset::Array{Int}
    starts::Array{Int}
    stops::Array{Int}
end

function Score(rows, cols; band=(lower=Inf,upper=Inf))
    num_elts = [(r == 0) ? 0 : min(r+band.upper, cols) - max(r-band.lower, 1) + 1 for r in 0:rows]
    offset   = round.(Int, cumsum(num_elts))
    starts   = round.([max(r-band.lower, 1) for r in 1:rows])
    stops    = round.([min(r+band.upper, cols) for r in 1:rows])
    return Score(rows, cols, band, zeros(offset[end]), offset, starts, stops)
end

# data stored in row-major form
# @inline index(S::Score, i, j) = (j-S.starts[i]+1) + S.offset[i]
@inline index(S::Score, i, j) = j + S.cols*(i-1)
@inline data(S::Score, i, j)  = S.data[index(S,i,j)]

# --------------------------------
# overloads of base operators

Base.size(S::Score)        = (S.rows, S.cols)
Base.ndims(T::Type{Score}) = 2

Base.BroadcastStyle(::Type{Score})              = Base.Broadcast.ArrayStyle{Score}()
Base.show(io::IO, S::Score)                     = Base.show(io, reshape(S.data, S.rows, S.cols))
Base.show(io::IO, ::MIME"text/plain", S::Score) = Base.show(io, reshape(S.data, S.rows, S.cols))
Base.display(S::Score)                          = Base.display(reshape(S.data, S.rows, S.cols))

function Base.lastindex(S::Score, d)
    if d == 1
        return S.rows
    elseif d == 2
        return S.cols
    else
        error("incorrect number of dimensions")
    end
end

# TODO: better (less-explicit) indexing
#       would be nice to keep as abstract ranges if we could, e.g. 1:2:10
#       right now we allocate memory every time
@inline Base.getindex(S::Score, inds...) = Base.getindex(S.data, index(S,r,c))
    # TODO: make this faster!
    # close(out)
    #       causes huge slow down
    # rows, cols = inds
    # ι = [index(S,r,c) for r in rows for c in cols]
    # if length(ι) == 1
    #     # if cols < S.starts[rows] || cols > S.stops[rows]
    #     #     println("out of bound access: (", rows, ",", cols, ")" )
    #     #     return -Inf
    #     # else
    #     return Base.getindex(S.data, ι[1])
    #     # end
    # else
    #     return Base.getindex(S.data, ι)
    # end
# end

function Base.setindex!(S::Score, X, inds...)
    i, j = inds
    Base.setindex!(S.data, X, index(S,i,j))
end

# ------------------------------------------------------------------------
# cigar functions

# alignment -> cigar
function cigar(seq₁::Array{UInt8}, seq₂::Array{UInt8})
    if length(seq₁) != length(seq₂)
        error("not an alignment")
    end

    aln = IOBuffer()
    M, I, D = 0, 0, 0
    for (c₁, c₂) in zip(seq₁, seq₂)
        @match (Char(c₁), Char(c₂)) begin
            ('-','-') => error("both columns are gaps")
            ('-', _ ) => begin
                if I > 0
                    write(aln, f"{I}I")
                    I = 0
                elseif M > 0
                    write(aln, f"{M}M")
                    M = 0
                end
                D += 1
            end
            ( _ ,'-') => begin
                if D > 0
                    write(aln, f"{D}D")
                    D = 0
                elseif M > 0
                    write(aln, f"{M}M")
                    M = 0
                end
                I += 1
            end
            ( _ , _ ) => begin
                if D > 0
                    write(aln, f"{D}D")
                    D = 0
                elseif I > 0
                    write(aln, f"{I}I")
                    I = 0
                end
                M += 1
            end
        end
    end

    if I > 0
        write(aln, f"{I}I")
        I = 0
    elseif M > 0
        write(aln, f"{M}M")
        M = 0
    elseif D > 0
        write(aln, f"{D}D")
        D = 0
    end

    return String(take!(aln))
end

function uncigar(cg::String)
    chan = Channel{Tuple{Int, Char}}(0)
    @async begin
        i₁, i₂ = 1, 1
        while i₁ <= length(cg)
            while isdigit(cg[i₂])
                i₂ += 1
            end
            put!(chan, (parse(Int,cg[i₁:i₂-1]),cg[i₂]))
            i₁ = i₂ + 1
            i₂ = i₁
        end
        close(chan)
    end

    return chan
end

# TODO: relax hardcoded cutoff
# TODO: relax hardcoded reliance on cigar suffixes. make symbols instead
# chunk alignment 
#
mutable struct Pos
    start::Int
    stop::Int
end

to_index(x::Pos) = x.start:x.stop
advance!(x::Pos) = x.start=x.stop

function homologous(alignment, qry::Array{UInt8}, ref::Array{UInt8}; maxgap=500)
    # ----------------------------
    # internal type needed for iteration
    
    SNPMap   = Dict{Int,UInt8}
    IndelMap = Dict{Int, Union{Array{UInt8, Int}}}
    
    qryₓ = Pos(1,1)
    refₓ = Pos(1,1)

    # ----------------------------
    # list of blocks and their mutations
    
    seq   = Array{UInt8}[]                                         # all blocks of alignment cigar
    pos   = NamedTuple{(:qry,:ref),Tuple(Maybe{Pos},Maybe{Pos})}[] # position corresponding to each block
    snp   = Union{SNPMap,Nothing}[]                                # snps of qry relative to ref
    indel = Union{IndelMap,Nothing}[]                              # indels of qry relative to ref

    # current block being constructed
    block = (
        len   = 0,
        seq   = IOBuffer(),
        snp   = SNPMap(),
        indel = IndelMap(),
    )

    # ----------------------------
    # list of blocks and their mutations
    
    function finalize_block!()
        if block.len <= 0
            @goto advance
        end
        push!(pos, (qry = copy(qryₓ), ref = copy(refₓ)))
        push!(seq, take!(block.seq))
        push!(snp, block.snp)
        push!(indel, block.indel)

        block.len   = 0
        block.snp   = SNPMap()
        block.indel = IndelMap()
        # block.seq is cleared by take! above
        
        @label advance
        advance!(qryₓ)
        advance!(refₓ)
    end

    # ----------------------------
    # main bulk of algorithm
    for (len, type) in alignment
        @match type begin
        # TODO: treat soft clips differently?
        'S' || 'H' => begin
            # TODO: implement
            error("need to implement soft/hard clipping")
        end
        'M' => begin
            x = Pos(refₓ.stop, refₓ.stop+len)
            y = Pos(qryₓ.stop, qryₓ.stop+len)

            write(block.seq, ref[x])

            for locus in findall(qry[y] .!= ref[x])
                block.snp[block.len+locus] = qry[qryₓ.stop+locus]
            end

            qryₓ.stop += len
            refₓ.stop += len
            block.len += len
        end
        'D' => begin
            if len >= maxgap
                finalize_block!()

                x = Pos(refₓ.start,refₓ.stop+len)

                push!(pos, (qry=nothing, ref=x))
                push!(seq, ref[x])
                push!(snp, nothing)
                push!(indel, nothing)

                advance!(refₓ)
            else
                block.indel[block.len] = len
                refₓ.stop += len
            end
        end
        'I' => begin
            if len >= maxgap
                finalize_block!()

                x = Pos(qryₓ.start,qryₓ.stop+len)

                push!(pos, (qry=x, ref=nothing))
                push!(seq, qry[x])
                push!(snp, nothing)
                push!(indel, nothing)

                advance!(qryₓ)
            else
                x = Pos(qryₓ.stop,qryₓ.stop+len)
                block.indel[block.len] = qry[x]
                qryₓ.stop += len
            end
        end
         _  => error("unrecognized cigar string suffix")
        end
    end

    finalize_block!()

    return seq, pos, snp, indel
end

# ------------------------------------------------------------------------
# alignment functions

# paf alignment pair
mutable struct Hit
    name::String
    length::Int
    start::Int
    stop::Int
    seq::Union{Array{UInt8},Nothing}
end

mutable struct Alignment
    qry::Hit
    ref::Hit
    matches::Int
    length::Int
    quality::Int
    orientation::Bool
    cigar::Union{String,Nothing}
    divergence::Union{Float64,Nothing}
    align::Union{Float64,Nothing}
end

# dynamic programming
function align(seq₁::Array{UInt8}, seq₂::Array{UInt8}, cost)
    if length(seq₁) > length(seq₂)
        seq₁, seq₂ = seq₂, seq₁
    end

    L₁, L₂ = length(seq₁)+1, length(seq₂)+1

    # initialize matrices
    score = (
         M = zeros(L₁,L₂), #Score(L₁,L₂,band=cost.band),
         I = zeros(L₁,L₂), #Score(L₁,L₂,band=cost.band),
         D = zeros(L₁,L₂), #Score(L₁,L₂,band=cost.band),
    )
    # score = (
    #      M = Score(L₁,L₂,band=cost.band),
    #      I = Score(L₁,L₂,band=cost.band),
    #      D = Score(L₁,L₂,band=cost.band),
    # )

    # NOTE: upper and lower could be flipped
    bound(x,d) = isinf(x) ? size(score.M,d) : min(size(score.M,d),x)
    for j in 1:bound(cost.band.upper,2)
        score.M[1,j] = cost.gap(j-1)
    end
    for i in 1:bound(cost.band.lower,1)
        score.M[i,1] = cost.gap(i-1)
    end

    score.I[1,1:bound(cost.band.upper,2)] .= -Inf
    score.D[1:bound(cost.band.lower,1),1] .= -Inf

    # fill in bulk
    for i in 2:L₁
        lb = round(Int, max(i-cost.band.lower, 2))
        ub = round(Int, min(i+cost.band.upper, L₂))
        for j in lb:ub
            score.I[i,j] = max(
               score.M[i-1,j]+cost.open,
               score.I[i-1,j]+cost.extend,
            )
            score.D[i,j] = max(
               score.M[i,j-1]+cost.open,
               score.D[i,j-1]+cost.extend,
            )
            score.M[i,j] = max(
               score.M[i-1,j-1]+cost.match(seq₁[i-1],seq₂[j-1]),
               score.I[i,j],
               score.D[i,j],
            )
        end
    end

    # backtrack
    a₁, a₂ = IOBuffer(), IOBuffer()
    i, j = size(score.M)
    while i > 1 && j > 1
        if score.M[i,j] == score.D[i,j]
            k = 1
            while j > k && ((score.M[i,j-k] + cost.gap(k)) != score.M[i,j])
                k += 1
            end
            write(a₁,'-'^k) 
            write(a₂,seq₂[j-1:-1:j-k]) 
            j -= k
        elseif score.M[i,j] == score.I[i,j]
            k = 1
            while i > k && ((score.M[i-k,j] + cost.gap(k)) != score.M[i,j])
                k += 1
            end
            write(a₁,seq₁[i-1:-1:i-k]) 
            write(a₂,"-"^k) 
            i -= k 
        else
            write(a₁,seq₁[i-1])
            write(a₂,seq₂[j-1])
            i -= 1
            j -= 1
        end
    end

    if i > 1
        write(a₁,seq₁[i-1:-1:1])
    end

    if j > 1
        write(a₂,seq₂[j-1:-1:1])
    end

    b₁ = reverse(take!(a₁))
    b₂ = reverse(take!(a₂))
    return b₁, b₂
end

include("static/watson-crick.jl")
function reverse_complement(seq::Array{UInt8})
    cmpl = [wcpair[nuc] for nuc in seq]
    reverse!(cmpl)
    return cmpl
end

cost = (
    open   = -6.0,
    extend = -1.0,
    band   = (
        lower = Inf,
        upper = Inf,
    ),
    gap    = k -> k == 0 ? 0 : cost.open + cost.extend*(k-1),
    match  = (c₁, c₂) -> 6.0*(c₁ == c₂) - 3.0,
)
# TODO: come up with a better function name
function enforce_cutoff!(a::Alignment, χ)
    δqₗ, δqᵣ = a.qry.start, a.qry.length - a.qry.stop
    δrₗ, δrᵣ = a.ref.start, a.ref.length - a.ref.stop

    s₁ = a.orientation ? a.qry.seq : reverse_complement(a.qry.seq)
    s₂ = a.ref.seq

    # left side of match
    if (0 < δqₗ <= χ) && (δrₗ == 0 || δrₗ > χ)
        a.qry.start = 0
        a.cigar     = string(δqₗ) * "I" * a.cigar
    elseif (δrₗ <= χ) && (δqₗ == 0 || δqₗ > χ)
        a.ref.start = 0
        a.cigar     = string(δrₗ) * "D" * a.cigar
    elseif (δqₗ <= χ) && (δrₗ <= χ)
        a₁, a₂ = align(s₁[1:δqₗ], s₂[1:δrₗ])
        cg     = cigar(a₁, a₂)

        a.qry.start = 0
        a.ref.start = 0
        a.cigar   = cg * a.cigar
        a.length += len(a₁)
    end

    # right side of match
    if (0 < δqᵣ <= χ) && (δrᵣ == 0 || δrᵣ > χ)
        a.qry.stop  = a.qry.length
        a.cigar     = a.cigar * string(δqₗ) * "I"
    elseif (δrᵣ <= χ) && (δqᵣ == 0 || δqᵣ > χ)
        a.ref.stop  = a.ref.length
        a.cigar     = a.cigar * string(δrₗ) * "D"
    elseif (δqᵣ <= χ) && (δrᵣ <= χ)
        a₁, a₂ = align(s₁[end-δqᵣ:end], s₂[end-δrᵣ:end])
        cg     = cigar(a₁, a₂)

        a.qry.start = a.qry.length
        a.ref.start = a.ref.length
        a.cigar   = a.cigar * cg
        a.length += len(a₁)
    end
end

# ------------------------------------------------------------------------
# io functions

log(msg) = println(stderr, msg)

# fasta sequence record
struct Record
    seq::Array{UInt8}
    name::String
    meta::String
end

name(r::Record) = isempty(r.meta) ? r.name : r.name * " " * r.meta 

NL = '\n'
Base.show(io::IO, rec::Record) = print(io, f">{rec.name} {rec.meta}{NL}{String(rec.seq[1:40])}...{String(rec.seq[end-40:end])}")

function read_fasta(io)
    chan = Channel{Record}(0)
    @async begin
        buf=IOBuffer()
        line=readline(io)
        while !isempty(line) && line[1] == '>'
            words      = split(line[2:end])
            name, meta = words[1], join(words[2:end], " ")

            line=readline(io)

            while !isempty(line) && line[1] != '>'
                write(buf,rstrip(line))
                line=readline(io)
            end
            put!(chan, Record(take!(buf), name, meta))
        end

        close(buf)
        close(chan)
    end

    return chan
end

function Base.show(io::IO, h::Hit)
    print(io, f"{h.name}[{h.length}]: ({h.start},{h.stop})")
end

function Base.show(io::IO, a::Alignment)
    print(io, f"qry: {a.qry}", '\t')
    print(io, f"ref: {a.ref}", '\t')
    print(io, f"polarity: {a.orientation}")
end

function read_paf(io)
    chan = Channel{Alignment}(0)

    int(x)   = parse(Int,x)
    float(x) = parse(Float64,x)
    last(x)  = split(x,':')[end]

    @async begin
        for row in eachline(io)
            elt = split(strip(row))
            cg = nothing
            dv = nothing
            as = nothing
            for x in elt[13:end]
                if startswith(x, "cg:")
                    cg = last(x)
                elseif startswith(x, "de:f")
                    dv = float(last(x))
                elseif startswith(x, "AS:i")
                    as = int(last(x))
                end
            end

            put!(chan, Alignment(Hit(elt[1],int(elt[2]),int(elt[3]),int(elt[4]),nothing),
                                 Hit(elt[6],int(elt[7]),int(elt[8]),int(elt[9]),nothing),
                                 int(elt[10]), int(elt[11]), int(elt[12]),
                                 elt[5] == "+",cg,dv,as))
        end
        close(chan)
    end

    return chan
end

# ------------------------------------------------------------------------
# string processing functions

function columns(s; nc=80)
    nr   = ceil(Int64, length(s)/nc)
    l(i) = 1+(nc*(i-1)) 
    r(i) = min(nc*i, length(s))
    rows = [String(s[l(i):r(i)]) for i in 1:nr]
    return join(rows,'\n')
end


function test()
    # println(">testing fasta parse...")
    # open("data/test.fna") do io
    #     for record in read_fasta(io)
    #         println(String(record.seq))
    #         println(record.seq)
    #     end
    # end
    # println(">done!")

    # println(">testing paf parse")
    # open("data/test.paf") do io
    #     for aln in read_paf(io)
    #         println(aln)
    #     end
    # end
    # println("done!")

    println(">testing cigar serialization...")
    s₁ = Vector{UInt8}("A-TCGT-GTCA-TAGC")
    s₂ = Vector{UInt8}("AGG-GTCGTCAGT-GC")
    cg = cigar(s₁, s₂)
    println("-->", cg)
    println(">done!")

    cost = (
        open   = -0.75,
        extend = -0.5,
        band   = (
            lower = Inf,
            upper = Inf,
        ),
        gap    = k -> k == 0 ? 0 : cost.open + cost.extend*(k-1),
        match  = (c₁, c₂) -> 2.0*(c₁ == c₂) - 1.0,
    )

    seq = (N) -> Vector{UInt8}(random_id(;len=N, alphabet=['A','C','G','T']))
    @benchmark align($seq(100), $seq(100), $cost)

    s = [ (seq(100), seq(100)) for i in 1:100 ]
    Profile.clear()
    @profile for (s₁, s₂) in s
        align(s₁, s₂, cost)
    end
    ProfileView.view()
    # println("1: ", String(a₁))
    # println("2: ", String(a₂))
end

end
