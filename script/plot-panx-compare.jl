using Statistics, StatsBase
using PanGraph
using Plots, ColorSchemes

include("plot-util.jl")

delta(X,x,len) = minimum(mod.(abs.(X .- x), len))

function compare(qry, ref; cutoff=150)
    data = Int[]
    len  = length(Graphs.sequence(qry))

    # sample empirical distribution
    for i in 1:length(qry.node)
        x, y  = qry.position[i], qry.position[i+1]
        depth = Graphs.depth(qry.node[i].block)
        if depth ≤ cutoff && length(qry.node[i]) ≥ 1000
            push!(data, delta(ref.position, x, len))
            push!(data, delta(ref.position, y, len))
        end
    end

    # generate null distribution (uniform sampling of fixed number of breakpoints)
    null = Array{Int,1}(undef,length(data))

    x = sort(sample(1:len, 1+length(data), replace=false))
    for i in 1:length(data)
        null[i] = delta(ref.position, x[i], len)
    end

    return (data, null)
end

function compare(graph)
    data = Int[]
    null = Int[]

    for (name, qry) in graph.pang.sequence
        d, n = compare(qry, graph.panx.sequence[name])
        append!(data,d); append!(null,n)
    end

    return data, null
end

function unmarshal(path)
    return (
        pang = open(Graphs.unmarshal,path.pang),
        panx = open(Graphs.unmarshal,path.panx),
    )
end

function main(path)
    data, null = path |> unmarshal |> compare

    p = plot(;
            xlabel="distance to nearest gene (bp)",
            ylabel="CDF",
            title="Klebsiella Pneumoniae",
            legend=:bottomright,
            xscale=:log10,
    )
    cdfplot!(data.+1;
        color=:black,
        linewidth=2,
        label="algorithm"
    )
    cdfplot!(null.+1;
        color=:red,
        linewidth=2,
        label="random null",
    )

    vline!([median(data)];
        label="",
        color=:black,
        linestyle=:dashdot,
    )
    vline!([median(null)];
        label="",
        color=:red,
        linestyle=:dashdot,
    )

    annotate!(20,  .85, "median≈70",  :black)
    annotate!(800, .25, "median≈300", :red)

    return p
end

if abspath(PROGRAM_FILE) == @__FILE__
    main((pang=ARGS[1], panx=ARGS[2]))
end
