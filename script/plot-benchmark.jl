using Statistics
using Plots, ColorSchemes

import Base: parse, +, /

struct Params
    N :: Int
    L :: Int
end

function parse(::Type{Params}, s::AbstractString)
    element = split(s,';')
    return Params(
        parse(Int,element[1]),
        parse(Int,element[2]),
    )
end

struct Duration
    hour   :: Int
    minute :: Int
    second :: Float64
end

function parse(::Type{Duration}, s::AbstractString)
    index = (
        o = 0,
        h = findfirst('h', s),
        m = findfirst('m', s),
        s = findfirst('s', s),
    )

    inc(::Nothing) = 1
    inc(x::Int64)  = x+1

    h = index.h === nothing ? 0 : parse(Int, s[inc(index.o):index.h-1])
    m = index.m === nothing ? 0 : parse(Int, s[inc(index.h):index.m-1])
    s = index.s === nothing ? 0 : parse(Float64, s[inc(index.m):index.s-1])

    return Duration(h,m,s)
end

seconds(t::Duration) = t.second + 60*t.minute + 3600*t.hour

function +(t1::Duration, t2::Duration)
    s = t1.second + t2.second
    m = t1.second + t2.second + (s÷60)
    h = t1.hour + t2.hour + (m÷60)

    return Duration(h, m%60, s%60)
end

function /(t::Duration, x::Int)
    h = t.hour ÷ x
    m = 60*(t.hour % x) + t.minute ÷ x
    s = 60*(t.minute % x) + t.second ÷ x

    m += s % 60
    s %= 60

    h += m % 60
    m %= 60

    return Duration(h, m, s)
end

function benchmarks(io)
    duration = Dict{Params,Array{Duration,1}}()
    for line in eachline(io)
        key, val = strip.(split(line, "=>"))
        param = parse(Params, key)
        dtime = parse(Duration, val)

        if param in keys(duration)
            push!(duration[param], dtime)
        else
            duration[param] = Duration[dtime]
        end
    end

    return duration
end

field(data, key) = Set(getfield.(keys(data), key)) |> collect |> sort
function plots(benchmark)
    L = field(benchmark, :L)
    N = field(benchmark, :N)

    colors = cgrad(:matter, length(N), categorical=true)
    p = plot(;
        legend = :topleft,
        xlabel = "bp/genome",
        ylabel = "time(sec)",
    )
    for (i,n) in enumerate(N)
        data = [seconds.(benchmark[Params(n,l)]) for l in L]
        μ, σ = mean.(data), std.(data)
        plot!(L, μ;
            ribbon=σ,
            fillalpha=0.25,
            color=colors[i],
            linewidth=2,
            xscale=:log10,
            yscale=:log10,
            label="population size=$(n)"
        )
    end

    plot!(L, 8e-4*L.*log10.(L), color=:black, linewidth=2, linestyle=:dashdot, label="log linear trend")

    p
end

function main(path, dest)
    plt = open(benchmarks, path) |> plots
    savefig(plt, dest)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(ARGS[1], ARGS[2])
end
