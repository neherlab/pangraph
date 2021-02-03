export Counter
export add!

Counter = Dict{String, Int}

function add!(c::Counter, s::String)
    if s ∈ keys(c)
        c[s] += 1
    else
        c[s]  = 1
    end
end
