using FileIO, JLD2
using Statistics, StatsBase
using PanGraph
using CairoMakie, ColorSchemes

# plot in headless mode -> allow for plotting on the server where no display is available
ENV["GKSwstype"] = "100"

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
    x = sort(sample(1:len, length(data)÷2, replace=false))
    null = Array{Int,1}(undef,length(x))

    for i in 1:length(x)
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
        "klebsiella_pneumoniae"  => "Klebsiella Pneumoniae",
        "mycobacterium_tuberculosis"  => "Mycobacterium Tuberculosis",
        "prochlorococcus_marinus" => "Prochlorococcus Mariunus",
        "escherichia_coli" => "Escherichia Coli",
        "helicobacter_pylori" => "Helicobacter Pylori",
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

    return pairs
end

function shorten(name)
    words = split(name)

    return "$(uppercase(words[1][1])). $(words[2])"
end

function main(paths)
    fig = Figure(font="Latin Modern Math", fontsize=26)
    axis = Axis(fig[1,1],
        xlabel="Species",
        ylabel="distance to nearest gene (bp)",
        xticks=(1:length(paths), map((p)->shorten(p.name), paths)),
        xticklabelrotation=π/6,
        yticks=(0:4, [L"10^0", L"10^1", L"10^2", L"10^3", L"10^4"]),
    )
    colors = cgrad(:Set1_5, length(paths), categorical=true)

    for (i,path) in enumerate(paths)
        plots!(axis, path, i, colors[i])
    end

    fig
end

function plots!(axis, path, i, color)
    exported = "$(dirname(path.pang))/export.jld2"
    data, null = path |> unmarshal |> compare

    x = vcat(fill(i,length(data)+length(null)))
    y = log10.(vcat(data, null) .+ 1)
    s = vcat(fill(:left, size(data)), fill(:right, size(null)))
    c = vcat(fill((color,1.0),size(data)), fill((color,0.5),size(null)))

    violin!(axis, x, y, side=s, color=c)
end


if abspath(PROGRAM_FILE) == @__FILE__
    plt = ARGS[2:end] |> group |> main
    CairoMakie.save("script/figs/$(ARGS[1]).png", plt, px_per_unit=2)
    CairoMakie.save("script/figs/$(ARGS[1]).pdf", plt)
end
