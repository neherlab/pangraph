
function log(msg...)
    println(stderr, msg...)
    flush(stderr)
end

function log_graph(G)
    msg = "### Graph ###"
    msg *= "\n--- Paths ---"
    # nodes = [n for (name, path) in G.sequence for n in path.node]
    for (name, path) in G.sequence
        msg *= "\n - $(name)"
        for n in path.node
            msg *= "\n\t$(n.block.uuid)|" * (n.strand ? "+" : "-")
            if n ∉ keys(G.block[n.block.uuid].mutate)
                msg *= " !!! node missing"
            end
        end
    end
    msg *= "\n--- Blocks ---"
    for (uuid, block) in G.block
        msg *= "\n - $(uuid) - length: $(length(block))"
        # for n in keys(block.mutate)
        #     if n ∉ nodes
        #         msg *= "\n\t ! node missing"
        #     end
        # end
    end
    log(msg)
end

Marginalize = Command(
    "marginalize",
    "pangraph marginalize <options> [arguments]",
    "computes all pairwise marginalizations of a multiple sequence alignment graph",
    """multiple sequence alignment accepted in formats: [json]""",
    [
        Arg(
            String,
            "output path",
            (short = "-o", long = "--output-path"),
            "path to directory where all pairwise marginalizations will be stored\n\tif empty, will skip this computation",
            "",
        ),
        Arg(
            Bool,
            "reduce paralog paths",
            (short = "-r", long = "--reduce-paralog"),
            "collapse coparallel paths through duplications",
            false,
        ),
        Arg(
            String,
            "isolates to project onto",
            (short = "-s", long = "--strains"),
            "collapse the graph to only blocks contained by paths of the given isolates.\n\tcomma seperated list, no spaces",
            "",
        ),
    ],
    (args) -> let
        path = parse(Marginalize, args)
        path = if (path === nothing || length(path) == 0)
            nothing
        elseif length(path) == 1
            path
        else
            return 2
        end

        graph = load(path, Marginalize)
        log("Loaded Graph")
        # log_graph(graph)
        names = collect(keys(graph.sequence))

        reduce = arg(Marginalize, "-r")
        output = arg(Marginalize, "-o")

        if length(output) > 0
            isdir(output) || mkpath(output)
            pairs = [(n₁, n₂) for n₁ in names for n₂ in names if n₁ < n₂]

            for (k, (name₁, name₂)) in enumerate(pairs)
                log("\r[$(k)/$(length(pairs))] $(name₁)-$(name₂)...")
                G = Graphs.copy(graph)
                # G = deepcopy(graph)
                # log_graph(G)
                log("- copy")
                Graphs.keeponly!(G, name₁, name₂)
                # log_graph(G)
                log("- keeponly!")
                if reduce
                    log("- reduce (start)")
                    changed = true
                    while changed
                        l = length(G.block)
                        Graphs.deparalog!(G)
                        Graphs.detransitive!(G)
                        changed = l != length(G.block)
                    end
                    log("- reduce (end)")
                else
                    log("- detransitive (start)")
                    # log_graph(G)
                    Graphs.detransitive!(G)
                    log("- detransitive (end)")
                end
                log("- finalize (start)")
                # recompute positions
                Graphs.finalize!(G)
                log("- finalize (end)")
                open("$(output)/$(name₁)-$(name₂).json", "w") do io
                    marshal(io, G; fmt = :json)
                end
            end
        end

        isolates = arg(Marginalize, "-s")
        if length(isolates) > 0
            names = split(isolates, ',')
            Graphs.keeponly!(graph, names...)

            if reduce
                changed = true
                while changed
                    l = length(graph.block)
                    Graphs.deparalog!(graph)
                    Graphs.detransitive!(graph)
                    changed = l != length(graph.block)
                end
            else
                Graphs.detransitive!(graph)
            end

            # recompute positions
            Graphs.finalize!(graph)
            marshal(stdout, graph; fmt = :json)
        end
    end,
)
