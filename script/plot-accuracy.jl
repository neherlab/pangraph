using Printf
using Plots, ColorSchemes
using JLD2, FileIO

include("plot-util.jl")

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

function plotcdfs(data, x, y; group="")
    i = size(data,1)÷2
    c = cgrad(:matter, size(data,2), categorical=true)
    p = plot(;
        legend = :bottomright,
        xscale = :log10,
        xlabel = "breakpoint misplacement (bp)",
        ylabel = "fraction of breakpoints",
        title  = length(group) > 0 ? "accuracy $(group)" : ""
    )

    for j in 1:size(data,2)
        cdfplot!(data[i,j].+1;
            linewidth = 2,
            color     = c[j],
            label     = @sprintf("%.1E", 5*y[j]),
         )
    end

    return p
end

const re = r"-([0-9]+).jld2"
function group(path)
    return match(re, path)[1]
end

function main(paths, destdir)
    data = [ load(path) for path in paths ]
    grid = [ plotgrid(entropy(datum)...;  group=group(path)) for (path,datum) in zip(paths,data) ]
    cdfs = [ plotcdfs(accuracy(datum)...; group=group(path)) for (path,datum) in zip(paths,data) ]

    save(plot, base, name) = savefig(plot, "$(destdir)/$(base)-$(name).png")
    base(name) = basename(name)[1:end-5]

    for (plt,path) in zip(grid,paths)
        save(plt,"heatmap",base(path))
    end

    for (plt,path) in zip(cdfs,paths)
        save(plt,"cdf",base(path))
    end

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
    main(ARGS)
end
