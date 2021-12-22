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

        data[key] !== nothing || continue
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

function plotgrid(data, x, y; group="")
    heatmap(string.(x), string.(10 .* y), data';
        xlabel = "HGT rate / genome / generation",
        ylabel = "pairwise diversity",
        title  = length(group) > 0 ? "intersection entropy $(group)" : ""
    )
end

function plotcdfs(data, x, y; group="")
    i = size(data,1)÷2
    c = cgrad(:matter, size(data,2), categorical=true)
    p = plot(;
        legend = :bottomright,
        xlabel = "avg. breakpoint misplacement / genome (bp)",
        ylabel = "fraction of breakpoints",
        title  = length(group) > 0 ? "accuracy $(group)" : ""
    )

    for j in 1:size(data,2)
        cdfplot!(data[i,j];
            linewidth = 2,
            color     = c[j],
            label     = "δ=$(10*y[j])",
         )
    end

    return p
end

const re = r"-([0-9]+).jld2"
function group(path)
    return match(re, path)[1]
end

function main(paths)
    data = [ load(path) for path in paths ]
    grid = [ plotgrid(entropy(datum)...;  group=group(path)) for (path,datum) in zip(paths,data) ]
    cost = [ plotgrid(entropy(datum)...;  group=group(path)) for (path,datum) in zip(paths,data) ]
    cdfs = [ plotcdfs(accuracy(datum)...; group=group(path)) for (path,datum) in zip(paths,data) ]

    return plot(grid...), plot(cdfs...)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(ARGS)
end
