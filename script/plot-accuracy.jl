using Printf
using LaTeXStrings
using Plots, ColorSchemes
using JLD2, FileIO

import CairoMakie

# plot in headless mode -> allow for plotting on the server where no display is available
ENV["GKSwstype"] = "100"

include("plot-util.jl")

α_mutrate_to_divergence = 12.8

Base.zero(x::Type{Array{Float64,1}}) = Float64[]

function unpack(key)
    elt = split(key,'/')
    return (
       hgt=parse(Float64,elt[1]),
       snp=parse(Float64,elt[2]),
       nit=parse(Int64,elt[3]),
       elt=elt[4],
    )
end

function params(data)
    hgt = Float64[]
    snp = Float64[]

    for key in keys(data)
        param = unpack(key)
        param.hgt ∈ hgt || push!(hgt, param.hgt)
        param.snp ∈ snp || push!(snp, param.snp)
    end

    sort!(hgt)
    sort!(snp)

    return hgt,snp
end


function collectentry(data, type, name, accumulate!, post!)
    hgt, snp = params(data)
    val = [ zero(type) for i in 1:length(hgt), j in 1:length(snp) ]
    num = zeros(Int,length(hgt),length(snp))
    for key in keys(data)
        param = unpack(key)
        param.elt == name || continue

        i = findfirst(hgt .== param.hgt)
        j = findfirst(snp .== param.snp)

        data[key] !== nothing   || continue
        !any(isnan.(data[key])) || continue

        accumulate!(val,i,j,data[key])
        num[i,j] += 1
    end
    return post!(val, num), hgt, snp
end

# specific getters
entropy(data)    = collectentry(data,Float64,"tiles",(arr,i,j,x)->arr[i,j]+=x, (arr,num)->arr./num)
accuracy(data)   = collectentry(data,Array{Float64,1},"costs",(arr,i,j,x)->append!(arr[i,j],x),(arr,num)->arr)
diversity(data)  = collectentry(data,Float64,"dists",(arr,i,j,x)->arr[i,j]+=x, (arr,num)->arr./num)
complexity(data) = collectentry(data,Float64,"nblks",(arr,i,j,x)->arr[i,j]+=x, (arr,num)->arr./num)

function plotgrid(data, x, y; group="", labels=true)
    heatmap(string.(x), [@sprintf("%.1E",5*Y) for Y in y], data';
        xlabel = labels ? "HGT rate / genome / generation" : "",
        ylabel = labels ? "pairwise diversity" : "",
        title  = length(group) > 0 ? "entropy $(group)" : "",
        legend = :none,
        clim   = (0,1.5),
    )
end

function plotcdfs(data, x, y; group="", fontsize=12, kwargs...)
    i = max(1,size(data,1)÷2)
    c = cgrad(:matter, size(data,2), categorical=true)
    p = plot(;
        legend = :bottomright,
        xscale = :log10,
        xlabel = "breakpoint misplacement (bp)",
        ylabel = "fraction of breakpoints",
        title  = length(group) > 0 ? "accuracy $(group)" : "",
        xtickfontsize  = fontsize,
        ytickfontsize  = fontsize,
        xguidefontsize = fontsize,
        yguidefontsize = fontsize,
        legendfontsize = fontsize,
        kwargs...
    )

    for j in 1:size(data,2)
        cdfplot!(data[i,j].+1;
            linewidth = 1,
            color     = c[j],
            label     = @sprintf("%.1E", α_mutrate_to_divergence*y[j]),
         )
    end

    return p
end

function publication(data, x, y)
    i = size(data,1)÷2
    c = cgrad(:matter, size(data,2), categorical=true)

    fig  = CairoMakie.Figure(font="Latin Modern Math", fontsize=26)
    axis = CairoMakie.Axis(fig[1,1],
                xscale=CairoMakie.log10,
                xlabel="breakpoint misplacement (bp)",
                ylabel="fraction of breakpoints",
                xticks=([1, 10, 100, 1000], [L"10^0", L"10^1", L"10^2", L"10^3"]),
                yticks=CairoMakie.LinearTicks(5)
    )
    CairoMakie.hidespines!(axis, :t, :r)

    for j in 1:size(data,2)
        points = sort(reduce(vcat, data[i,j] for i in 1:size(data,1))).+1
        CairoMakie.lines!(axis, points , range(0,1,length(points)),
            color = c[j],
            label = @sprintf("%.1e", α_mutrate_to_divergence*y[j])
        )
    end

    CairoMakie.axislegend("Avg. pairwise diversity", position = :rb, nbanks=3)

    return fig
end

const re = r"/accuracy-([^/]+)\.jld2"
function group(path)
    return match(re, path)[1]
end

function main(path, destdir)
    data = load(path)
    grid = plotgrid(entropy(data)...;  group=group(path))
    cdfs = plotcdfs(accuracy(data)...; group=group(path))

    save(base, name) = savefig("$(destdir)/$(base)-$(name).png")
    save(plot, base, name) = savefig(plot, "$(destdir)/$(base)-$(name).png")

    base(name) = basename(name)[1:end-5]

    save(grid,"heatmap",base(path))
    save(cdfs,"cdf",base(path))

    # publication plot
    fig = publication(accuracy(data)...)
    CairoMakie.save("$(destdir)/paper-$(base(path)).png", fig, px_per_unit=2)
#=
    TODO: return one figure with correct layout
    l = @layout[Plots.grid(1,2) a{0.05w}]
    return (
        plot(grid...,
            heatmap((0:0.01:1).*ones(101,1);
                    legend=:none,
                    xticks=:none,
                    yticks=(1:10:101, [@sprintf("%1.1f",1.5*x) for x in 0:0.1:1]),
             ); layout=l
        ),
        plot(cdfs...),
        plot(grid[1], cdfs[1])
    )
=#
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(ARGS[1], ARGS[2])
end
