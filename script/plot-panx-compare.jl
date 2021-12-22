using Plots, ColorSchemes
using PanGraph

include("plot-util.jl")

function compare(qry, ref)
end

function main(path)
    graph = (
        pang = open(Graphs.unmarshal,path.pang),
        panx = open(Graphs.unmarshal,path.panx),
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    main((pang=ARGS[1], panx=ARGS[2]))
end
