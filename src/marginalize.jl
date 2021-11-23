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
        "path to direcotry where output files will be stored",
        "./tmp"
    ),
    Arg(
        Bool,
        "reduce paralog paths",
        (short="-r", long="--reduce-paralog"),
        "collapse coparallel paths through duplications",
        false,
    ),
   ],
   (args) -> let
       path = parse(Marginalize, args)
       path === nothing && return 2
       length(path) > 1 && return 2

       graph  = load(path, Marginalize)
       names  = collect(keys(graph.sequence))
       output = arg(Marginalize, "-o")

       isdir(output) || mkpath(output)

       pairs = [ (n₁,n₂) for n₁ in names for n₂ in names if n₁ < n₂ ]

       reduce = arg(Marginalize, "-r")
       Threads.@threads for (name₁, name₂) in pairs
           G = deepcopy(graph)
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
           Graphs.finalize!(G)
           open("$(output)/$(name₁)-$(name₂).json", "w") do io
               marshal(io, G; fmt=:json)
           end
       end
   end
)
