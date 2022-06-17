##
using Printf
using LaTeXStrings
using ColorSchemes
using JLD2, FileIO
using Plots
import CairoMakie
##

ENV["GKSwstype"] = "100"

Base.zero(x::Type{Array{Float64,1}}) = Float64[]

function unpack(key)
    elt = split(key, '/')
    return (
        hgt = parse(Float64, elt[1]),
        snp = parse(Float64, elt[2]),
        nit = parse(Int64, elt[3]),
        elt = elt[4],
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

    return hgt, snp
end


function collectentry(data, type, name, accumulate!, post!)
    hgt, snp = params(data)
    val = [zero(type) for i = 1:length(hgt), j = 1:length(snp)]
    num = zeros(Int, length(hgt), length(snp))
    for key in keys(data)
        param = unpack(key)
        param.elt == name || continue

        i = findfirst(hgt .== param.hgt)
        j = findfirst(snp .== param.snp)

        data[key] !== nothing || continue
        !any(isnan.(data[key])) || continue

        accumulate!(val, i, j, data[key])
        num[i, j] += 1
    end
    return post!(val, num), hgt, snp
end

# specific getters
diversity(data) = collectentry(
    data,
    Float64,
    "dists",
    (arr, i, j, x) -> arr[i, j] += x,
    (arr, num) -> arr ./ num,
)

mut_dens(data) = collectentry(
    data,
    Float64,
    "mut_dens",
    (arr, i, j, x) -> arr[i, j] += x,
    (arr, num) -> arr ./ num,
)

const rex = r"/accuracy-([^/]+)\.jld2"
function group(path)
    return match(rex, path)[1]
end

function linreg(fit_arr::Vector{Tuple{Float64,Float64}})::Float64
    x, y = [[f[i] for f in fit_arr] for i = 1:2]
    α = sum(x .* y) / sum(x .^ 2)
    return α
end

##

function comparison_plot(paths)
    ms = Dict("minimap10" => :xcross, "minimap20" => :cross, "mmseqs" => :circle)

    p = plot(
        xlabel = "mutation rate",
        ylabel = "avg. block divergence (inferred graph)",
        legend = :topleft,
    )
    fit_arr = Tuple{Float64,Float64}[]
    for (n, path) in enumerate(paths)
        data = load(path)
        path_id = group(path)

        path_label = Dict(
            "mmseqs" => "mmseqs2",
            "minimap10" => "minimap2 (asm10)",
            "minimap20" => "minimap2 (asm20)",
        )

        D, h, s = diversity(data)
        I, J = size(D)

        colors = cgrad(:thermal, I, categorical = true)

        for i = 1:I
            label = i == 1 ? path_label[path_id] : false
            scatter!(
                s,
                D[i, :],
                label = label,
                marker = ms[path_id],
                markercolor = colors[i],
            )

            label = n == 3 ? "hgt = $(h[i])" : false
            plot!(s, D[i, :], label = label, color = colors[i], linestyle = :dot)
            for j = 1:J
                s[j] <= 0.002 && push!(fit_arr, (s[j], D[i, j]))
            end
        end
    end

    α = linreg(fit_arr)
    s = 0:0.001:0.006
    plot!(
        s,
        s .* α,
        color = :gray,
        linestyle = :dash,
        label = "lin. reg. (α = $(round(α, digits=1)))",
    )

    return p
end


function main(figname, paths)
    sort!(paths)

    p = comparison_plot(paths)
    savefig(figname)

end

##

# paths = [
#     "script/synthetic_data/results/accuracy-$kernel.jld2" for
#     kernel in ["minimap10", "minimap20", "mmseqs"]
# ]
# main("script/figs/paper-accuracy-comparison.pdf", paths)


##

if abspath(PROGRAM_FILE) == @__FILE__
    main(ARGS[1], ARGS[2:end])
end

##
