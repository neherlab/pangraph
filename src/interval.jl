module Intervals

import Base:
    +, -, \,
    <, >, ==, !=, ≤,
    ∈, ∩, ∪, ⊆, ⊇, ~,
    max, min,
    isless, isequal, isdisjoint,
    iterate, isempty, length, copy, getindex, to_index,
    show

export Interval, IntervalSet
export containing

# -------------------------------------------------------------------------
# atomic interval

# represents an atomic (closed, open) interval
struct Interval{T<:Real}
    lo :: T
    hi :: T
    Interval{T}(lo,hi) where T <: Real = hi <= lo ? throw(ArgumentError("invalid interval hi=$hi ≤ low=$lo")) : new(lo, hi)
end

# ---------------------------------
# constructors

Interval(lo::T, hi::T) where T <: Real = Interval{T}(lo,hi)
Interval(x::T) where T <: Real = Interval{T}(x,x)
Interval(t::Tuple) = Interval(t...)

copy(a::Interval) = Interval(a.lo, a.hi)

# ---------------------------------
# repl properties

show(io::IO, a::Interval) = print(io, "[$(a.lo),$(a.hi))")

# ---------------------------------
# operations

# collection/

length(a::Interval) = a.hi - a.lo
to_index(a::Interval{T}) where T <: Integer = a.lo:(a.hi-1)

# logical

==(a::Interval, b::Interval) = a.lo == b.lo && a.hi == b.hi
≤(a::Interval, b::Interval)  = a.lo ≤ b.lo && a.hi ≤ b.hi
<(a::Interval, b::Interval)  = a.lo < b.lo && a.hi < b.hi

isless(a::Interval, b::Interval)  = a < b
isequal(a::Interval, b::Interval) = a == b

# set relationships
∈(x::T, a::Interval{T}) where T <: Real = a.lo ≤ x ≤ a.hi

⊆(a::Interval, b::Interval) = b.lo ≤ a.lo && a.hi ≤ b.hi
⊂(a::Interval, b::Interval) = (a==b) ? false : a ⊆ b
⊇(a::Interval, b::Interval) = b ⊆ a
⊃(a::Interval, b::Interval) = b ⊂ a

# arithmetic
+(a::Interval) = a
-(a::Interval) = Interval(-a.lo, -a.hi)

+(a::Interval{T}, x::T) where T <: Real = Interval{T}(a.lo+x, a.hi+x)
-(a::Interval{T}, x::T) where T <: Real = Interval{T}(a.lo-x, a.hi-x)

+(x::T, a::Interval{T}) where T <: Real = +a+x
-(x::T, a::Interval{T}) where T <: Real = -a+x

isdisjoint(a::Interval, b::Interval) = b.hi < a.lo || a.hi < b.lo

function ∩(a::Interval, b::Interval)
    try
        return Interval(max(a.lo,b.lo), min(a.hi,b.hi))
    catch err
        if isa(err, ArgumentError)
            return nothing # intervals were disjoint
        end

        throw(err) # unexpected error
    end
end

function \(a::Interval, b::Interval)
    ab = a ∩ b

    isnothing(ab) && return [a]                         # disjoint
    ab == a       && return typeof(a)[]                 # a ⊆ b
    a.lo == ab.lo && return [Interval(ab.hi, a.hi)]     # a ≥ b
    a.hi == ab.hi && return [Interval(a.lo, ab.lo)]     # a ≤ b

    return [Interval(a.lo, b.lo), Interval(b.hi, a.hi)] # a ⊃ b
end

# -------------------------------------------------------------------------
# interval sets

# ---------------------------------
# helpers

# NOTE: assumes (closed, open) intervals
function merge!(I::Array{Interval{T}}) where T <: Real
    i = 1
    while i < length(I)
        if I[i].hi < I[i+1].lo
            i += 1
            continue
        end

        I[i] = Interval(I[i].lo, I[i+1].hi)
        deleteat!(I, i+1)
    end
end

# ---------------------------------
# typedef

struct IntervalSet{T<:Real} 
    Is  :: Array{Interval{T},1}
    min :: T
    max :: T

    function IntervalSet{T}(Is, min, max) where T <: Real
        # ensures atomic intervals are disjoint
        𝕀 = Interval{T}[]
        for I in Is
            n = length(𝕀)
            W = [I]
            for i ∈ 1:n
                W = [x for w in W for x in w \ 𝕀[i]]
            end

            append!(𝕀, W)
        end

        sort!(𝕀)
        merge!(𝕀)

        return new(𝕀, min, max)
    end
end

# ---------------------------------
# constructors

IntervalSet(min::T, max::T, itr)                    where T <: Real    = IntervalSet{T}(itr, min, max)
IntervalSet(min::T, max::T, Is::Array{Interval{T}}) where T <: Real    = IntervalSet{T}(Is, min, max)
IntervalSet(min::T, max::T, I::Interval{T})         where T <: Real    = IntervalSet{T}([I], min, max)

IntervalSet(min::T, max::T)                         where T <: Real    = IntervalSet{T}(Interval{T}[], min, max)
IntervalSet(Is::Array{Interval{T}})                 where T <: Real    = IntervalSet{T}(Is, typemin(T), typemax(T))
IntervalSet(I::Interval{T})                         where T <: Real    = IntervalSet{T}([I], typemin(T), typemax(T))

IntervalSet(min, max, Is::Tuple{T,T}...)            where T <: Real    = IntervalSet{T}((Interval(I) for I in Is), min, max)
IntervalSet(Is::Tuple{T,T}...)                      where T <: Real    = IntervalSet{T}((Interval(I) for I in Is), typemin(T), typemax(T))

copy(I::IntervalSet{T}) where T <: Real = IntervalSet{T}(copy(I.Is), I.min, I.max)

# ---------------------------------
# repl properties

show(io::IO, Is::IntervalSet) = print(io, Is.Is)

# ---------------------------------
# collection properties

length(I::IntervalSet)  = length(I.Is)

iterate(I::IntervalSet)        = length(I) > 0 ? (I.Is[1], 2) : nothing
iterate(I::IntervalSet, state) = (state <= length(I)) ? (I.Is[state], state+1) : nothing

getindex(I::IntervalSet, i) = getindex(I.Is, i)
getindex(A::T, I::IntervalSet) where T <: AbstractArray = [a for i in I for a in A[i]]

min(I::IntervalSet) = I.Is[1].lo
max(I::IntervalSet) = I.Is[end].hi

# ---------------------------------
# operations

# logical
isempty(I::IntervalSet) = length(I) == 0

# arithmetic

function ~(I::IntervalSet) 
    starts = (i == 0 ? I.min : I[i].hi for i in 0:length(I))
    stops  = ((i<=length(I)) ? I[i].lo : I.max for i in 1:(length(I)+1))

    return IntervalSet(I.min, I.max, (Interval(start, stop) for (start, stop) in zip(starts,stops) if start < stop))
end

\(I::IntervalSet{T}, j::Interval{T})    where T <: Real = IntervalSet{T}((x for i ∈ I for x in i \ j), I.min, I.max)
\(I::IntervalSet{T}, J::IntervalSet{T}) where T <: Real = foldl(\, J; init=I)

# NOTE: extra generator voodoo avoids creating a list with nothings and then filtering them out 
∩(I::IntervalSet{T}, j::Interval{T}) where T <: Real = 
    IntervalSet{T}(
        ( x for i ∈ I for x in [i ∩ j] if !isnothing(x) ),
        I.min,
        I.max
    )
∩(j::Interval, I::IntervalSet) = I ∩ j

∩(I::IntervalSet{T}, J::IntervalSet{T}) where T <: Real = 
    IntervalSet{T}(
       sort([ x for i ∈ I for j ∈ J for x in [i ∩ j] if !isnothing(x) ]),
       min(I.min,J.min),
       max(I.max,J.max)
    )

function ∪(I::IntervalSet, J::IntervalSet)
    IJ = copy(I) 
    K  = J \ I
    for k ∈ K
        ι = searchsortedfirst(IJ.Is, k)
        insert!(IJ.Is, ι, k)
        merge!(IJ.Is)
    end

    return IJ
end

∪(I::IntervalSet, j::Interval) = I ∪ IntervalSet(j) 

isdisjoint(A::IntervalSet, b::Interval) = all(isdisjoint(a,b) for a in A)

function containing(A::IntervalSet{T}, x::T) where T <: Real
    for a in A
        if x ∈ a
            return a
        end
    end

    return nothing
end

function containing(A::IntervalSet, b::Interval)
    for a in A
        if b ⊆ a
            return a
        end
    end

    return nothing
end

# -------------------------------------------------------------------------
# tests

function test()
    i = Interval(1, 10)
    j = Interval(0, 11)

    @show i
    @show j

    @show i \ j
    @show j \ i

    I = IntervalSet(0, 50, (0, 4), (2, 5), (7, 8), (10,15), (17,19))
    J = IntervalSet(0, 50, (3, 4), (6, 8), (11,14), (16,18))

    @show I
    @show J

    @show I \ i
    @show I \ j

    @show I \ J

    @show I ∩ j

    @show I ∪ J

    @show ~I

    return nothing
end

end
