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

function group(args)
    i = 1
    sorted = sort(args)
    species = Dict(
        "kleb"  => "Klebsiella Pneumoniae",
        "myco"  => "Mycobacterium Tuberculosis",
        "proc"  => "Prochlorococcus Mariunus",
        "ecoli" => "Escherichia Coli",
        "helio" => "Helicobacter Pylori",
    )

    pairs = NamedTuple{(:pang,:panx,:name),Tuple{String,String,String}}[]
    while i ≤ length(sorted)
        push!(pairs, (
            pang=sorted[i],
            panx=sorted[i+1],
            name=species[basename(dirname(sorted[i]))],
        ))
        i += 2
    end

    @show pairs
    return pairs
end

function main(paths)
    figure = plot(;
            xlabel="distance to nearest gene (bp)",
            ylabel="CDF",
            legend=:bottomright,
            xscale=:log10,
    )
    colors = cgrad(:matter, length(paths), categorical=true)

    for (i,path) in enumerate(paths)
        plots!(path, colors[i])
    end

    figure
end

function plots!(path, color)
    data, null = path |> unmarshal |> compare

    cdfplot!(data.+1;
        color=color,
        linewidth=2,
        label="$(path.name) data"
    )
    cdfplot!(null.+1;
        color=color,
        linewidth=2,
        linestyle=:dashdot,
        label="$(path.name) null",
    )

    #=
    vline!([median(data)];
        label="",
        color=color,
        linestyle=:dashdot,
    )
    vline!([median(null)];
        label="",
        color=color,
        linestyle=:dashdot,
    )

    annotate!(20,  .85, "median≈70",  :black)
    annotate!(800, .25, "median≈300", :red)
    =#
end

function save(plt)
    savefig(plt, "figs/panx-compare.png")
end

if abspath(PROGRAM_FILE) == @__FILE__
    ARGS |> group |> main |> save
end
