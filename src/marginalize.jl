Marginalize = Command(
   "marginalize",
   "pangraph marginalize <options> [arguments]",
   "computes all pairwise marginalizations of a multiple sequence alignment graph",
   """multiple sequence alignment accepted in formats: [json]""",
   [
    Arg(
        String,
        "output path",
        (short="-o", long="--output-path"),
        "path to directory where all pairwise marginalizations will be stored\n\tif empty, will skip this computation",
        ""
    ),
    Arg(
        Bool,
        "reduce paralog paths",
        (short="-r", long="--reduce-paralog"),
        "collapse coparallel paths through duplications",
        false,
    ),
    Arg(
        String,
        "isolates to project onto",
        (short="-s", long="--strains"),
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

       graph  = load(path, Marginalize)
       names  = collect(keys(graph.sequence))

       reduce = arg(Marginalize, "-r")
       output = arg(Marginalize, "-o")

       if length(output) > 0
           isdir(output) || mkpath(output)
           pairs  = [ (n₁,n₂) for n₁ in names for n₂ in names if n₁ < n₂ ]

           Threads.@threads for (name₁, name₂) in pairs
               G = Graphs.copy(graph)
               Graphs.keeponly!(G, name₁, name₂)
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
               open("$(output)/$(name₁)-$(name₂).json", "w") do io
                   marshal(io, G; fmt=:json)
               end
           end
       end

       isolates = arg(Marginalize, "-s")
       if length(isolates) > 0
           names = split(isolates,',')
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
           marshal(stdout, graph; fmt=:json)
       end
   end
)
