module Mash

using ..Graphs: sequence

const map = UInt64[
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
]

const maxU64 = typemax(UInt64)

struct Minimizer
    value    :: UInt64
    position :: UInt64
end

# transliteration of the invertible hash function found in minimap
function hash(x::UInt64, mask::UInt64)::UInt64
    x = (~x + (x << 21)) & mask
    x = x ⊻ x >> 24
    x = (x + (x << 3) + (x << 8)) & mask
    x = x ⊻ x >> 14
    x = (x + (x << 2) + (x << 4)) & mask
    x = x ⊻ x >> 28
    x = (x + (x << 31)) & mask
    return x
end

function sketch(seq::Array{UInt8}, k::Int, w::Int, id::Int)
    (k < 0 || k > 32)  && error("k='$(k)' must be ∈ [0,32]")
    (w < 0 || w > 255) && error("w='$(w)' must be ∈ [0,255]")

    fwd :: UInt64 = 0
    rev :: UInt64 = 0

    mask  :: UInt64 = (1 << (2*k)) - 1;
    shift :: UInt64 = 2*(k-1);

    min = Minimizer(maxU64, maxU64)
    minimizer = Minimizer[]
    window = fill(Minimizer(maxU64,maxU64), w)

    l = 0
    bi, mi = 1, 1
    for (locus,nuc) in enumerate(seq)
        c = map[nuc+1] # need to offset as Julia is 1-indexed
        new = if c ≥ 4
            l = 0
            Minimizer(maxU64,maxU64)
        else
            fwd = ((fwd<<2) | c) & mask
            rev = ((rev>>2) | ((3⊻c) << shift))
            l += 1

            if l ≥ k
                pos = (UInt64(id) << 32) | (UInt64(locus) << 1)
                if fwd ≤ rev
                    Minimizer(hash(fwd,mask),pos)
                else
                    Minimizer(hash(rev,mask),pos|1)
                end
            else
                Minimizer(maxU64,maxU64)
            end
        end

        window[bi] = new
        if (l == w+k-1) && (min.value != maxU64)
            for i in (bi+1):w
                if min.value == window[i].value && min.position != window[i].position
                    push!(minimizer, window[i])
                end
            end
            for i in 1:bi
                if min.value == window[i].value && min.position != window[i].position
                    push!(minimizer, window[i])
                end
            end
        end

        if new.value < min.value
            if l ≥ w+k && min.value != maxU64
                push!(minimizer, min)
            end
            min = new
            mi = bi
        elseif bi == mi
            if (l ≥ w+k-1) && min.value != maxU64
                push!(minimizer, min)
            end

            min = Minimizer(maxU64, min.position)
            for i in (bi+1):w
                if window[i].value < min.value
                    mi  = i
                    min = window[i]
                end
            end

            for i in 1:bi
                if window[i].value < min.value
                    mi  = i
                    min = window[i]
                end
            end

            if (l ≥ w+k-1) && min.value != maxU64
                for i in (bi+1):w
                    if min.value == window[i].value && min.position != window[i].position
                        push!(minimizer, window[i])
                    end
                end
                for i in 1:bi
                    if min.value == window[i].value && min.position != window[i].position
                        push!(minimizer, window[i])
                    end
                end
            end
        end

        bi += 1
        if bi > w
            bi = 1
        end
    end

    if min.value != maxU64
        push!(minimizer, min)
    end

    return minimizer
end

function distance(graphs...; k=20, w=200)
    sequences = Dict(seq for graph in graphs for seq in sequence(graph))

    names = collect(keys(sequences))
    seqs  = [sequences[name] for name in names]

    minimizers = Minimizer[]
    for (i,seq) in enumerate(seqs)
        append!(minimizers, sketch(Array{UInt8}(seq), k, w, i))
    end

    @time sort!(minimizers, by=(m)->m.value)
    l = 1
    while l ≤ length(minimizers)
        l = r
        while r ≤ length(minimizers) && minimizers[r].value == minimizers[l].value
            r += 1
        end
    end

    exit(1)
end

end
