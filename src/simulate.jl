module Simulation

using Random, Distributions

# ------------------------------------------------------------------------
# types

struct Params
    pop_size :: Int
    seq_len  :: Int

    snp_rate :: Float64
    hgt_rate :: Float64
    del_rate :: Float64
    inv_rate :: Float64
end
Params(;pop_size=100, seq_len=Int(1e6), snp_rate=1e-6, hgt_rate=0, del_rate=0, inv_rate=0) = Params(pop_size, seq_len, snp_rate, hgt_rate, del_rate, inv_rate)

# bitpacked: 30 bytes(ancestor) | 30 bytes (location) | 3 bytes (mutation) | 1 byte strand
Sequence = Array{UInt64,1}
# states of mutuation
# 0** -> no mutation
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
    s[from:to] .⊻= 1
    reverse!(s, from, to)
end

function model(param::Params)
    parent = Array{Int}(undef, param.pop_size)
    indel = Normal(param.seq_len, param.seq_len/10)

    return function(population::Array{Sequence}, offspring::Array{Sequence})
        rand!(parent, 1:param.pop_size)

        for n in 1:param.pop_size
            copy!(offspring[n], population[parent[n]])

            # test to mutate
            L = length(offspring[n])
            μ = rand(Poisson(L*param.snp_rate))
            if μ > 0
                for locus in rand(1:L, μ)
                    mutate!(offspring[n], locus)
                end
            end

            # test to pick up donation
            if rand() ≤ param.hgt_rate
                donor = population[rand(1:param.pop_size)]

                from = rand(1:length(donor))
                area = Int(rand(indel))
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
            if rand() ≤ param.del_rate
                from = rand(1:length(offspring[n]))
                area = Int(rand(indel))
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
            if rand() ≤ param.inv_rate
                from = rand(1:length(offspring[n]))
                area = Int(rand(indel))
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

function test()
    N = 50
    L = Int(1e6)
    evolve! = model(Params(;pop_size=N,seq_len=L,snp_rate=5e-6))

    isolate   = population(fill(L, N))
    offspring = [ Array{UInt64,1}(undef,L) for _ in 1:N ]

    for _ in 1:2*N
        evolve!(isolate, offspring)
        isolate, offspring = offspring, isolate
    end

    for n in 1:N
        @show Int(mode(unpack.ancestor.(isolate[n])))
        @show Int(sum(unpack.snp.(isolate[n])))
    end

    return isolate
end

end
