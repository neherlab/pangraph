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


# dynamic programming
# NOTE: assumes not banded for now
#       will need to have a matrix-like object but with bands
function align(seq₁::Array{Char}, seq₂::Array{Char}, cost)
    L₁, L₂ = length(seq₁)+1, length(seq₂)+1

    # initialize matrices
    score = (
         M = zeros(L₁, L₂),
         I = zeros(L₁, L₂),
         D = zeros(L₁, L₂),
    )

    for j in 1:size(score.M,2)
        score.M[1,j] = cost.gap(j-1)
    end
    for i in 1:size(score.M,1)
        score.M[i,1] = cost.gap(i-1)
    end

    score.I[1,:] .= -Inf
    score.D[:,1] .= -Inf

    # fill in bulk
    for i in 2:L₁
        for j in 2:L₂
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

    s₁ = collect("MASSivvvvvvey")
    s₂ = collect("MISSivey")

    cost = (
        open   = -1,
        extend = -0.5,
        gap = k -> k == 0 ? 0 : cost.open + cost.extend*(k-1),
        match = (c₁, c₂) -> 2.0*(c₁ == c₂) - 1.0,
    )

    a₁, a₂ = align(s₁, s₂, cost)
    println("1: ", String(a₁))
    println("2: ", String(a₂))
end

end
