push!(LOAD_PATH, "../src/")
using Documenter, PanGraph

makedocs(
    sitename = "PanGraph.jl",
    authors  = "Nicholas Noll, Marco Molari, Richard Neher",
    modules  = [PanGraph],
    pages    = [
        "Home"    => "index.md",
        "Library" => [
            "lib/align.md",
            "lib/block.md",
            "lib/graph.md",
            "lib/mash.md",
            "lib/minimap.md",
            "lib/node.md",
            "lib/path.md",
            "lib/simulate.md",
        ],
        "Command Line" => [
            "cli/build.md",
            "cli/generate.md",
        ]
    ]
)
