module Utility

using FStrings
using StatsBase

import Base.Threads.@spawn

export random_id, log
export read_fasta, name
export read_paf

# random string of fixed length
function random_id(;length=10)
    alphabet = ['A','B','C','D','E','F','G','H','I','J','K','L','M',
                'N','O','P','Q','R','S','T','U','V','W','X','Y','Z']
    return join(sample(alphabet, length))
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

struct Hit
    name::String
    length::Int
    start::Int
    stop::Int
end

function Base.show(io::IO, h::Hit)
    print(io, f"{h.name}[{h.length}]: ({h.start},{h.stop})")
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

function columns(s; nc=80)
    nr   = ceil(Int64, length(s)/nc)
    l(i) = 1+(nc*(i-1)) 
    r(i) = min(nc*i, length(s))
    rows = [String(s[l(i):r(i)]) for i in 1:nr]
    return join(rows, '\n')
end

function test()
    println("testing fasta parse...")
    open("data/test.fna") do io
        for record in read_fasta(io)
            println(String(record.seq))
            println(columns(record.seq))
        end
    end
    println("done!")

    println("testing paf parse")
    open("data/test.paf") do io
        for aln in read_paf(io)
            println(aln)
        end
    end
    println("done!")
end

end
