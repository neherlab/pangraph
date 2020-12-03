module Utility

using FStrings, Match
using StatsBase
using Infiltrator

import Base.Threads.@spawn

export random_id, log
export read_fasta, name
export read_paf

# ------------------------------------------------------------------------
# random functions

# random string of fixed length
function random_id(;length=10)
    alphabet = ['A','B','C','D','E','F','G','H','I','J','K','L','M',
                'N','O','P','Q','R','S','T','U','V','W','X','Y','Z']
    return join(sample(alphabet, length))
end

# ------------------------------------------------------------------------
# cigar/alignment functions

# paf alignment pair
struct Hit
    name::String
    length::Int
    start::Int
    stop::Int
end

struct Alignment
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
function index(S::Score, i, j) 
    return (j-S.starts[i]+1) + S.offset[i]
end
data(S::Score, i, j)  = S.data[index(S,i,j)]

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
function Base.getindex(S::Score, inds...)
    rows, cols = inds
    ι = [index(S,r,c) for r in rows for c in cols]
    if length(ι) == 1
        if cols < S.starts[rows] || cols > S.stops[rows]
            println("out of bound access: (", rows, ",", cols, ")" )
            return -Inf
        else
            return Base.getindex(S.data, ι[1])
        end
    else
        return Base.getindex(S.data, ι)
    end
end

function Base.setindex!(S::Score, X, inds...)
    i, j = inds
    Base.setindex!(S.data, X, index(S,i,j))
end

# dynamic programming
function align(seq₁::Array{Char}, seq₂::Array{Char}, cost)
    if length(seq₁) > length(seq₂)
        seq₁, seq₂ = seq₂, seq₁
    end

    L₁, L₂ = length(seq₁)+1, length(seq₂)+1

    # initialize matrices
    score = (
         M = Score(L₁,L₂,band=cost.band),
         I = Score(L₁,L₂,band=cost.band),
         D = Score(L₁,L₂,band=cost.band),
    )

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
            println("I")
            score.I[i,j] = max(
               score.M[i-1,j]+cost.open,
               score.I[i-1,j]+cost.extend,
            )
            println("D")
            score.D[i,j] = max(
               score.M[i,j-1]+cost.open,
               score.D[i,j-1]+cost.extend,
            )
            println("M")
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

    return collect(reverse(String(take!(a₁)))), collect(reverse(String(take!(a₂))))
end

function cigar(seq₁::Array{Char}, seq₂::Array{Char})
    if length(seq₁) != length(seq₂)
        error("not an alignment")
    end

    aln = IOBuffer()
    M, I, D = 0, 0, 0
    for (c₁, c₂) in zip(seq₁, seq₂)
        @match (c₁, c₂) begin
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

include("static/watson-crick.jl")
function reverse_complement(seq::Array{Char})
    cmpl = [wcpair[nuc] for nuc in seq]
    reverse!(cmpl)
    return cmpl
end

function extend!(a::Alignment, χ)
    δqₗ, δqᵣ = a.qry.start, a.qry.length - a.qry.stop
    δrₗ, δrᵣ = a.ref.start, a.ref.length - a.ref.stop

    # left side of match
    if     (δqₗ <= χ) && (δrₗ > χ || δrₗ == 0)
        a.qry.start = 0
        a.cigar     = δqₗ * "I" * a.cigar
    elseif (δrₗ <= χ) && (δqₗ > χ || δqₗ == 0)
        a.ref.start = 0
        a.cigar     = δrₗ * "D" * a.cigar
    elseif (δqₗ <= χ) && (δrₗ <= χ)
    end

    # right side of match
end

# ------------------------------------------------------------------------
# io functions

log(msg) = println(stderr, msg)

# fasta sequence record
struct Record
    seq::Array{Char}
    name::String
    meta::String
end

name(r::Record) = isempty(r.meta) ? r.name : r.name * " " * r.meta 

NL = '\n'
Base.show(io::IO, rec::Record) = print(io, f">{rec.name} {rec.meta}{NL}{String(rec.seq[1:40])}...{String(rec.seq[end-40:end])}")

function read_fasta(io)
    chan = Channel{Record}(0)
    @spawn begin
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
    last(x)  = split(x)[end]

    @spawn begin
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

            put!(chan, Alignment(Hit(elt[1],int(elt[2]),int(elt[3]),int(elt[4])),
                                 Hit(elt[6],int(elt[7]),int(elt[8]),int(elt[9])),
                                 int(elt[10]), int(elt[11]), int(elt[12]),
                                 elt[5] == "+",cg,dv,as))
        end

        close(chan)
    end

    return chan
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
    s₁ = collect("A-TCGT-GTCA-TAGC")
    s₂ = collect("AGG-GTCGTCAGT-GC")
    cg = cigar(s₁, s₂)
    println("-->", cg)
    println(">done!")

    s₁ = collect("MSS")
    s₂ = collect("MIIISS")
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

    a₁, a₂ = align(s₁, s₂, cost)
    println("1: ", String(a₁))
    println("2: ", String(a₂))
end

end
