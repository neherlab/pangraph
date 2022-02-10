push!(LOAD_PATH, "../src/")
using Documenter, PanGraph

makedocs(
    sitename = "PanGraph.jl",
    authors = "Nicholas Noll, Marco Molari, Richard Neher",
    modules = [PanGraph],
    pages = [
        "Home" => "index.md",
        "Tutorial" => [
            "tutorials/tutorial_1.md",
            "tutorials/tutorial_2.md",
            "tutorials/tutorial_3.md",
            "tutorials/tutorial_4.md",
        ],
        "Library" => [
            "lib/align.md",
            "lib/block.md",
            "lib/edge.md",
            "lib/graph.md",
            "lib/gfa.md",
            "lib/mash.md",
            "lib/minimap.md",
            "lib/node.md",
            "lib/path.md",
            "lib/simulate.md",
            "lib/utility.md",
        ],
        "Command Line" => [
            "cli/build.md",
            "cli/export.md",
            "cli/generate.md",
            "cli/marginalize.md",
            "cli/polish.md",
        ],
    ],
)
