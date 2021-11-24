module Simulation

using Statistics, Random, Distributions

include("interval.jl")
using .Intervals

include("graph.jl")
using .Graphs

# ------------------------------------------------------------------------
# types

struct Rates
    snp :: Float64
    hgt :: Float64
    del :: Float64
    inv :: Float64
end
Rates(;snp=0, hgt=0, del=0, inv=0) = Rates(snp, hgt, del, inv)

struct Params
    N    :: Int
    L    :: Int
    σₗ   :: Int
    rate :: Rates
end

Params(; N=100, L=Int(1e6), σₗ=Int(1e5), snp=0, hgt=0, del=0, inv=0) = Params(N, L, σₗ, Rates(snp, hgt, del, inv))

# bitpacked: 30 bytes(ancestor) | 30 bytes (location) | 3 bytes (mutation) | 1 byte strand
Sequence = Array{UInt64,1}
# states of mutuation
# 0** -> no mutation, ** ignored
# 1** -> mutation given by **
const shift = (
    ancestor = 34,
    location = 4,
    mutation = 1,
    strand   = 0,
)

mask(n::Int) = (64 ≤ n || n ≤ 0) ? 0 : UInt64((1 << n)-1)

const unpack = (
    ancestor = (x::UInt64) -> (x>>shift.ancestor) & mask(30),
    location = (x::UInt64) -> (x>>shift.location) & mask(30),
    mutation = (x::UInt64) -> (x>>shift.mutation) & mask(2),
    strand   = (x::UInt64) -> x & mask(1),

    # helpers
    snp = (x::UInt64) -> (x >> 3) & 1
)

ancestor(n::Int,L::Int)   = [(UInt64(n) << shift.ancestor) | (UInt64(l) << shift.location) for l in 1:L]
population(L::Array{Int}) = [ ancestor(n,len) for (n,len) in enumerate(L)]

# ------------------------------------------------------------------------
# evolution operators

mutate!(s::Sequence, at::Int) = s[at] |= ((rand(0:3) << 1) | 8)
delete!(s::Sequence, from::Int, to::Int) = splice!(s, from:to, [])
insert!(acceptor::Sequence, donor::Sequence, at::Int) = splice!(acceptor, at:at+1, donor)
function invert!(s::Sequence, from::Int, to::Int)
    s[from:to] .⊻= 1        # complement
    reverse!(s, from, to)   # reverse
end

function model(param::Params)
    parent = Array{Int}(undef, param.N)
    indel = Normal(param.L, param.σₗ)

    int(x) = Int(round(x))

    return function(population::Array{Sequence}, offspring::Array{Sequence})
        rand!(parent, 1:param.N)

        for n in 1:param.N
            copy!(offspring[n], population[parent[n]])

            # test to mutate
            L = length(offspring[n])
            μ = rand(Poisson(L*param.rate.snp))
            if μ > 0
                for locus in rand(1:L, μ)
                    mutate!(offspring[n], locus)
                end
            end

            # test to pick up donation
            if rand() ≤ param.rate.hgt
                donor = population[rand(1:param.N)]

                from = rand(1:length(donor))
                area = int(rand(indel))
                if area > length(offspring[n])
                    d  = area - length(offspring[n])
                    to = (from + d - 1) % length(donor) + 1

                    if to > from
                        donor = donor[from:to]
                    else
                        donor = vcat(donor[from:length(donor)], donor[1:to])
                    end

                    site = rand(1:length(offspring[n]))
                    insert!(offspring[n], donor, site)
                end
            end

            # test to delete
            if rand() ≤ param.rate.del
                from = rand(1:length(offspring[n]))
                area = int(rand(indel))
                if area < length(offspring[n])
                    d  = length(offspring[n]) - area
                    to = (from + d - 1) % length(offspring[n]) + 1

                    if to > from
                        delete!(offspring[n], from, to)
                    else
                        delete!(offspring[n], from, length(offspring[n]))
                        delete!(offspring[n], 1, to)
                    end
                end
            end

            # test to invert
            if rand() ≤ param.rate.inv
                from = rand(1:length(offspring[n]))
                area = int(rand(indel))
                d  = abs(area - length(offspring[n]))
                to = (from + d - 1) % length(offspring[n]) + 1

                if to > from
                    invert!(offspring[n], from, to)
                else
                    invert!(offspring[n], from, length(offspring[n]))
                    invert!(offspring[n], 1, to)
                end
            end
        end
    end
end

function pancontig!(s::Sequence, ancestor::Dict{Int,Array{Interval}})
	isolate = NamedTuple{(:loci, :strand, :ancestor), Tuple{Interval,Bool,Int}}[]

	# helper functions
	function push_ancestor!(a, loci)
		if a ∈ keys(ancestor)
			atoms = ancestor[a]
			disjoint = IntervalSet(loci) \ IntervalSet(atoms[1].lo, atoms[end].hi, atoms)
			overlaps = [ x for atom in ancestor[a] for x in partition(atom,loci) ]
			if length(disjoint) == 0
				ancestor[a] = overlaps
			else
				ancestor[a] = sort(vcat(overlaps, disjoint.Is))
			end
        else
			ancestor[a] = [loci]
        end
	end

	function push_isolate!(a, loci, strand)
		push!(isolate, (
			loci=loci,
			strand=strand,
			ancestor=a
		))
	end

	function interval(left, right)
		left < right && return Interval(left, right+1)
		return Interval(right, left+1)
	end

	function increments(prev, curr, strand)
		strand && return prev+1 == curr # forward strand
		return prev-1 == curr # reverse strand
	end

	l  = 1
    id = unpack.ancestor(s[1])
    start = last = unpack.location(s[1])
	strand = unpack.strand(s[1]) == 0

	for (r,state) in enumerate(s[2:end])
        locus = unpack.location(state)
        ident = unpack.ancestor(state)
		polar = unpack.strand(state) == 0

        # still in tile
		if (polar == strand) && (ident == id) && increments(last, locus, strand)
            last = locus
            continue
        end

		push_ancestor!(id, interval(start, last))
		push_isolate!(id, interval(start, last), strand)

		l  = r+1
        id = ident
        start = last = locus
		strand = polar
    end

	# hanging block
	if start != last
		push_ancestor!(id, interval(start, last))
		push_isolate!(id, interval(start, last), strand)
	end

	return isolate
end

function pancontigs(isolates::Array{Sequence})
	ancestors = Dict{Int, Array{Interval}}()
	isolates  = [ pancontig!(isolate, ancestors) for isolate in isolates ]

	return isolates, ancestors
end

function nodes(isolate, ancestors, blocks)
	ancestor = ancestors[isolate.ancestor]
	subpath  = Graphs.Node[]
	for interval in ancestor
		if isolate.loci ⊇ interval
			block = blocks[(isolate.ancestor,interval)]
			node  = Graphs.Node(block, isolate.strand)

			block.mutate[node] = Graphs.SNPMap()
			block.insert[node] = Graphs.InsMap()
			block.delete[node] = Graphs.DelMap()

			push!(subpath, node)
		end
	end
	return subpath
end

const ntab = UInt8['A', 'C', 'G', 'T' ]
function nucleotide(sequence::Array{Sequence}, ancestor::Array{Array{UInt8,1},1})
	function transform(x::UInt64)
		if unpack.snp(x) == 1
			return ntab[unpack.mutation(x)+1]
		else
			ident  = unpack.ancestor(x)
			locus  = unpack.location(x)
			return ancestor[ident][locus]
		end
	end
	transform(s::Sequence) = transform.(s)

	return transform.(sequence)
end

function run(evolve!, time::Int, initial::Array{Array{UInt8,1},1}; graph=false)
	sequence  = population(length.(initial))
	offspring = [ Array{UInt64,1}(undef,maximum(length.(initial))) for _ in 1:length(initial) ]

    for _ in 1:time
        evolve!(sequence, offspring)
        sequence, offspring = offspring, sequence
    end

	if graph
		isolates, ancestors = pancontigs(sequence)

		# build the evolved pangraph from ancestral intervals
		blocks = Dict(
			(ancestor,interval) => Graphs.Block(initial[ancestor][interval])
			for (ancestor,intervals) in ancestors for interval in intervals
		)
		paths = [
			let
				node = [ n for interval in isolate for n in nodes(interval, ancestors, blocks) ]
				Graphs.Path("isolate_$(i)", node, 0, true, [])
			end for (i,isolate) in enumerate(isolates)
		]

		G = Graphs.Graph(Dict(pair(b) for b in values(blocks)), Dict(pair(p) for p in paths))
		Graphs.detransitive!(G)
		Graphs.finalize!(G)

		return nucleotide(sequence,initial), (isolates, ancestors, G)
	end

	return nucleotide(sequence,initial), nothing
end

randseq(len::Int) = rand(UInt8.(['A','C','G','T']), len)

function test()
    N = 10
    L = Int(1e4)
	evolve! = model(Params(;N=N,L=L,snp=1e-3,hgt=0,inv=0))
	ancestors = [randseq(L) for _ in 1:N]
	return run(evolve!, 20, ancestors; graph=true)
end


end
