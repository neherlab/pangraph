module Simulation

using Random, Distributions

include("src/interval.jl")
using .Intervals

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
    mutation = (x::UInt64) -> (x>>shift.mutation) & mask(3),
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

                    continue
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

                    continue
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

                continue
            end
        end
    end
end

function pancontig!(s::Sequence, lengths::Int, ancestor::Dict{Int,IntervalSet})
    start = last = unpack.location(s[0])
    id = unpack.ancestor(s[0])
    for state in s[2:end]
        locus = unpack.location(state)
        ident = unpack.ancestor(state)

        # still in tile
        if ((locus == last + 1 || locus == 1 && last == lengths[indent])
        &&  (ident == id))
            last = locus
            continue
        end

        if id ∈ ancestor
            ancestor[id] = add(ancestor[id], Interval(start, locus))
        else
            ancestor[id] = IntervalSet((start,locus))
        end
        start = last = locus
        id = ident
    end
end

function test()
    N = 50
    L = Int(1e6)
    r = 5e-1
    μ = 0
    evolve! = model(Params(;N=N,L=L,snp=μ,hgt=r))

    isolate   = population(fill(L, N))
    offspring = [ Array{UInt64,1}(undef,L) for _ in 1:N ]

    for _ in 1:N
        evolve!(isolate, offspring)
        isolate, offspring = offspring, isolate
    end

    for n in 1:N
        @show length(isolate[n])
        @show unique(Int.(unpack.ancestor.(isolate[n])))
        @show Int(sum(unpack.snp.(isolate[n])))
    end

    return isolate
end

end
