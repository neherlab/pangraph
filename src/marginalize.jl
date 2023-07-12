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
        Arg(
            Bool,
            "check that genomes are preserved",
            (short = "-t", long = "--test"),
            "consistency check that verifies that genomes in the output graphs are identical\n\tto genomes in the input graphs.",
            false,
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
        names = collect(keys(graph.sequence))

        reduce = arg(Marginalize, "-r")
        output = arg(Marginalize, "-o")
        test_flag = arg(Marginalize, "-t")

        function verify(G₀, G₋)
            # verifies that after the marginalization all genomes from
            # G₋ are exactly equal to their initial value in in G₀
            for (name, path) in G₋.sequence
                seq₀ = sequence(G₀.sequence[name], shift = true)
                seq₋ = sequence(path, shift = true)
                @assert seq₀ == seq₋
            end
        end

        function log(msg...)
            println(stderr, msg...)
            flush(stderr)
        end

        if length(output) > 0
            isdir(output) || mkpath(output)
            pairs = [(n₁, n₂) for n₁ in names for n₂ in names if n₁ < n₂]

            Threads.@threads for (name₁, name₂) in pairs
                Γ = Graphs.copy_graph(graph)
                Graphs.keeponly!(Γ, name₁, name₂)
                if reduce
                    changed = true
                    while changed
                        l = length(Γ.block)
                        Graphs.deparalog!(Γ)
                        Graphs.detransitive!(Γ)
                        changed = l != length(Γ.block)
                    end
                else
                    Graphs.detransitive!(Γ)
                end

                # recompute positions
                Graphs.finalize!(Γ)

                test_flag && verify(graph, Γ)

                open("$(output)/$(name₁)-$(name₂).json", "w") do io
                    marshal(io, Γ; fmt = :json)
                end
            end
        end

        isolates = arg(Marginalize, "-s")
        if length(isolates) > 0

            G = Graphs.copy_graph(graph)
            names = split(isolates, ',')
            Graphs.keeponly!(G, names...)

            if reduce
                changed = true
                while changed
                    l = length(G.block)
                    Graphs.deparalog!(G)
                    Graphs.detransitive!(G)
                    changed = l != length(G.block)
                end
            else
                Graphs.detransitive!(G)
            end

            # recompute positions
            Graphs.finalize!(G)

            test_flag && verify(graph, G)

            marshal(stdout, G; fmt = :json)
        end
    end,
)
