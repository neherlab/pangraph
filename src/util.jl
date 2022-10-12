module Utility

using Rematch
using StatsBase

# internal modules
using ..Intervals, ..Nodes

import ..Graphs:
    reverse_complement, reverse_complement!,
    SNPMap, InsMap, DelMap

import ...PanGraph: Alignment, Hit

# exports
export random_id, log
export cigar, uncigar
export hamming_align

export contiguous_trues
export make_consensus, alignment_alleles

export write_fasta, read_fasta, name
export read_paf
export lock_semaphore

# ------------------------------------------------------------------------
# multithreading resource allocation

function lock_semaphore(f::Function, s::Base.Semaphore)
    Base.acquire(s)
    try
        return f()
    finally
        Base.release(s)
    end
end

# ------------------------------------------------------------------------
# random functions


# random string of fixed length
"""
	random_id(;len=10, alphabet=UInt8[])

Generate a random string of length `len` drawn from letters in `alphabet`.
"""
function random_id(;len=10, alphabet=UInt8[])
    if length(alphabet) == 0
        alphabet = ['A','B','C','D','E','F','G','H','I','J','K','L','M',
                    'N','O','P','Q','R','S','T','U','V','W','X','Y','Z']
    end
    return join(sample(alphabet, len))
end

# helper functions w/ common functionality
"""
	write_fasta(io::IO, name, seq)

Output a single FASTA record with sequence `seq` and name `name` to IO stream `io`.
"""
function write_fasta(io::IO, name, seq)
    write(io::IO, '>', name, '\n')
    write(io::IO, columns(seq), '\n')
end

# ------------------------------------------------------------------------
# banded alignment data structure

"""
	struct Score <: AbstractArray{Float64,2}
		rows::Int
		cols::Int
		band::NamedTuple{(:lower, :upper)}
		data::Array{Float64}
		offset::Array{Int}
		starts::Array{Int}
		stops::Array{Int}
	end

Store information about a banded pairwise alignment.
"""
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
"""
	cigar(seq₁::Array{UInt8}, seq₂::Array{UInt8})

Given two sequences, `seq₁` and `seq₂`, perform a pairwise banded alignment and return the cigar string of alignment.
"""
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
                    write(aln, "$(I)I")
                    I = 0
                elseif M > 0
                    write(aln, "$(M)M")
                    M = 0
                end
                D += 1
            end
            ( _ ,'-') => begin
                if D > 0
                    write(aln, "$(D)D")
                    D = 0
                elseif M > 0
                    write(aln, "$(M)M")
                    M = 0
                end
                I += 1
            end
            ( _ , _ ) => begin
                if D > 0
                    write(aln, "$(D)D")
                    D = 0
                elseif I > 0
                    write(aln, "$(I)I")
                    I = 0
                end
                M += 1
            end
        end
    end

    if I > 0
        write(aln, "$(I)I")
        I = 0
    elseif M > 0
        write(aln, "$(M)M")
        M = 0
    elseif D > 0
        write(aln, "$(D)D")
        D = 0
    end

    return String(take!(aln))
end

uncigar(cg::T) where T <: AbstractArray{Tuple{Int,Char}} = cg

"""
	uncigar(cg::String)

Return an interator over intervals of alignment defined by cigar string `cg`.
"""
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

# ------------------------------------------------------------------------
# alignment functions

# dynamic programming
"""
	align(seq₁::Array{UInt8}, seq₂::Array{UInt8}, cost::Score)

Perform a pairwise alignment using Needleman-Wunsch style dynamic programming between
`seq₁` and `seq₂` given `cost`. The cost is defined by the Score structure.
"""
function align(seq₁::Array{UInt8}, seq₂::Array{UInt8}, cost)
    flip = if length(seq₁) > length(seq₂)
        seq₁, seq₂ = seq₂, seq₁
        true
    else
        false
    end

    L₁, L₂ = length(seq₁)+1, length(seq₂)+1

    # initialize matrices
    score = (
         M = zeros(L₁,L₂), #Score(L₁,L₂,band=cost.band),
         I = zeros(L₁,L₂), #Score(L₁,L₂,band=cost.band),
         D = zeros(L₁,L₂), #Score(L₁,L₂,band=cost.band),
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
            write(a₂,'-'^k)
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
        write(a₂,'-'^(i-1))
    elseif j > 1
        write(a₂,seq₂[j-1:-1:1])
        write(a₁,'-'^(j-1))
    elseif j > 1 && i > 1
        error("invalid backtraversal")
    end

    b₁ = reverse(take!(a₁))
    b₂ = reverse(take!(a₂))
    return flip ? (b₂, b₁) : (b₁, b₂)
end

"""
	hamming_align(qry::Array{UInt8,1}, ref::Array{UInt8,1})

Perform a simple alignment of `qry` to `ref` by minimizing hamming distance.
Useful for fast, approximate alignments of small sequences.
"""
function hamming_align(qry::Array{UInt8,1}, ref::Array{UInt8,1})
    matches = [
        let
            len = min(length(ref)-x+1, length(qry))
            sum(qry[1:len] .== ref[x:x+len-1])
        end
    for x ∈ 1:length(ref) ]

    return argmax(matches)
end

include("static/watson-crick.jl")
# XXX: This is the major allocator for us.
#      Will require serious thought as to how to not create so much garbage
"""
	reverse_complement(seq::Array{UInt8})

Return a newly allocated sequence array that is the reverse complement of `seq`.
"""
reverse_complement(seq::Array{UInt8}) = UInt8[wcpair[nuc+1] for nuc in reverse(seq)]

"""
	reverse_complement!(hit::Hit)

Reverse complement the qry of Hit in place.
"""
function reverse_complement!(hit::Hit)
    if hit.seq !== nothing
        hit.seq = reverse_complement(hit.seq)
    end

    start, stop = hit.start, hit.stop

    hit.start = 1 + (hit.length - stop)
    hit.stop  = hit.length - (start - 1)
end

function reverse_complement!(aln::Alignment)
    reverse_complement!(aln.qry)
    aln.orientation = ~aln.orientation
end

# ------------------------------------------------------------------------
# alignment modification

"""
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

`cost` are the default dynamic alignment parameters used.
"""
const cost = (
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
"""
	enforce_cutoff!(a::Alignment, χ)

Ensure that the alignment `a` does not have insertion or deletion segments larger than `χ`.
Return the list of segments created by parsing the alignment such that all segments are larger than `χ`.
"""
function enforce_cutoff!(a::Alignment, χ)
    δqₗ, δqᵣ = a.qry.start-1, a.qry.length - a.qry.stop
    δrₗ, δrᵣ = a.ref.start-1, a.ref.length - a.ref.stop

    s₁ = a.orientation ? a.qry.seq : reverse_complement(a.qry.seq)
    s₂ = a.ref.seq

    # left side of match
    if (0 < δqₗ ≤ χ) && (δrₗ == 0 || δrₗ > χ)
        a.qry.start = 1
        pushfirst!(a.cigar, (δqₗ, 'I'))
    elseif (0 < δrₗ ≤ χ) && (δqₗ == 0 || δqₗ > χ)
        a.ref.start = 1
        pushfirst!(a.cigar, (δrₗ, 'D'))
    elseif (0 < δqₗ ≤ χ) && (0 < δrₗ <= χ)
        a₁, a₂ = align(s₁[1:δqₗ], s₂[1:δrₗ], cost)
        cg     = collect(uncigar(cigar(a₁, a₂)))

        a.qry.start = 1
        a.ref.start = 1

        prepend!(a.cigar, cg)
        a.length += length(a₁)
    end

    # right side of match
    if (0 < δqᵣ ≤ χ) && (δrᵣ == 0 || δrᵣ > χ)
        a.qry.stop  = a.qry.length
        push!(a.cigar, (δqᵣ, 'I'))
    elseif (0 < δrᵣ ≤ χ) && (δqᵣ == 0 || δqᵣ > χ)
        a.ref.stop  = a.ref.length
        push!(a.cigar, (δrᵣ, 'D'))
    elseif (0 < δqᵣ ≤ χ) && (δrᵣ <= χ)
        a₁, a₂ = align(s₁[end-δqᵣ+1:end], s₂[end-δrᵣ+1:end], cost)
        cg     = collect(uncigar(cigar(a₁, a₂)))

        a.qry.stop = a.qry.length
        a.ref.stop = a.ref.length

        append!(a.cigar, cg)
        a.length += length(a₁)
    end
end

# ------------------------------------------------------------------------
# insertion -> intervals

# ------------------------------------------------------------------------
# io functions

log(msg...) = println(stderr, msg...)

# fasta sequence record
"""
	struct Record
		seq::Array{UInt8}
		name::String
		meta::String
	end

A record obtained when parsing a single entry of a FASTA file.
"""
struct Record
    seq  :: Array{UInt8}
    name :: String
    meta :: String
end

name(r::Record) = isempty(r.meta) ? r.name : r.name * " " * r.meta

NL = '\n'
Base.show(io::IO, rec::Record) = print(io, ">$(rec.name) $(rec.meta)$(NL)$(String(rec.seq[1:40]))...$(String(rec.seq[end-40:end]))")

"""
	read_fasta(io::IO)

Parse a FASTA file from IO stream `io`.
Return an iterator over all records.
"""
function read_fasta(io::IO)
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
    print(io, "$(h.name)[$(h.length)]: ($(h.start),$(h.stop))")
end

function Base.show(io::IO, a::Alignment)
    print(io, "qry: $(a.qry)", '\t')
    print(io, "ref: $(a.ref)", '\t')
    print(io, "polarity: $(a.orientation)")
end

"""
	read_paf(io::IO)

Parse a PAF file from IO stream `io`.
Return an iterator over all pairwise alignments.
"""
function read_paf(io::IO)
    int(x)   = parse(Int,x)
    float(x) = parse(Float64,x)
    last(x)  = split(x,':')[end]

    return [ let
        elt = split(strip(row))
        cg = nothing
        dv = nothing
        as = nothing
        for x in elt[13:end]
            if startswith(x, "cg:")
                cg = last(x)
            elseif startswith(x, "de:f")
                dv = last(x) |> float
            elseif startswith(x, "gi:f")
                dv = last(x) |> float
                dv = 1. - (dv/100.)
            elseif startswith(x, "AS:i")
                as = last(x) |> int
            end
        end

        # XXX: important: shift PAF entry to be 1 indexed and right-inclusive to match julia
        Alignment(
            Hit(
                elt[1],
                int(elt[2]),   # length
                int(elt[3])+1, # start
                int(elt[4]),   # stop
                nothing,
            ),
            Hit(
                elt[6],
                int(elt[7]),   # length
                int(elt[8])+1, # start
                int(elt[9]),   # stop
                nothing,
            ),
            int(elt[10]),
            int(elt[11]),
            int(elt[12]),
            elt[5] == "+",
            String(cg),
            dv,
            as
        )
        end for row in eachline(io) ]
end

"""
    read_mmseqs2(io::IO)

Parse a simil-PAF file produced by mmseq2 from IO stream `io`.
Return an iterator over all pairwise alignments.
"""
function read_mmseqs2(io::IO)
    int(x)   = parse(Int,x)
    float(x) = parse(Float64,x)
    last(x)  = split(x,':')[end]

    return [ let
        elt = split(strip(row))

        # XXX: important: shift PAF entry to be 1 indexed and right-inclusive to match julia
        q1 = int(elt[3])
        q2 = int(elt[4])
        qstart, qend, strand = q1 < q2 ? (q1, q2, true) : (q2, q1, false)

        Alignment(
            Hit(
                elt[1],
                int(elt[2]),   # length
                qstart, # start
                qend,   # stop
                nothing,
            ),
            Hit(
                elt[6],
                int(elt[7]),   # length
                int(elt[8]), # start
                int(elt[9]),   # stop
                nothing,
            ),
            int(elt[10]),
            int(elt[11]),
            int(elt[12]),
            strand,
            String(elt[13]),
            1. - float(elt[14]),
            int(elt[15])
        )
        end for row in eachline(io) ]
end

# ------------------------------------------------------------------------
# string processing functions

"""
	columns(s; nc=80)

Partition string `s` into an array of strings such that no string is longer than `nc` characters.
"""
function columns(s; nc=80)
    L = length(s)
    nr   = ceil(Int64, L/nc)
    l(i) = 1+(nc*(i-1)) 
    r(i) = min(nc*i, L)
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

    s = [ (seq(100), seq(100)) for i in 1:100 ]
    # println("1: ", String(a₁))
    # println("2: ", String(a₂))
end

function contiguous_trues(x)
    intervals = Interval{Int}[]

    l      = 1
    inside = false

    for (r,σ) ∈ enumerate(x)
        if σ
            if !inside
                l = r
                inside = true
            end

            continue
        end

        if inside
            inside = false
            push!(intervals, Interval(l, r))
        end
    end

    # check for case that extends to right boundary
    if inside
        push!(intervals, Interval(l, length(x)+1))
    end

    return IntervalSet(intervals)
end

# ------------------------------------------------------------------------
# alignment processing functions

make_consensus(alignment) = [mode(view(alignment,i,:)) for i in 1:size(alignment,1)]

function alignment_alleles(ref, aln, nodes)
    isdiff = aln .!= ref
    refdel = ref .== UInt8('-')
    alndel = aln .== UInt8('-')

    δ = (
        snp =   isdiff .& .~refdel .& .~alndel,
        del = .~refdel .&   alndel,
        ins =   refdel .& .~alndel,
    )

    coord   = cumsum(.!refdel)
    refgaps = contiguous_trues(refdel)
    gaps    = Dict{Int,Int}(coord[gap.lo] => length(gap) for gap in refgaps)

    mutate = Dict{Node,SNPMap}(
            node => SNPMap(
                   coord[l] => aln[l,i]
                for l in findall(δ.snp[:,i])
            )
        for (i,node) in enumerate(nodes)
    )

    delete = Dict{Node,DelMap}(
            node => DelMap(
                      del.lo => length(del)
                for del in contiguous_trues(δ.del[.~refdel,i])
             )
        for (i,node) in enumerate(nodes)
    )

    Δ(I) = (R = containing(refgaps, I)) == nothing ? 0 : I.lo - R.lo
    insert = Dict{Node,InsMap}(
            node => InsMap(
                      (coord[ins.lo],Δ(ins)) => aln[ins,i]
                for ins in contiguous_trues(δ.ins[:,i])
             )
        for (i,node) in enumerate(nodes)
    )

    sequence = ref[.~refdel]

    return (
        gaps     = gaps,
        mutate   = mutate,
        delete   = delete,
        insert   = insert,
        sequence = sequence,
    )
end

end
