Polish = Command(
   "polish",
   "pangraph polish <options> [pangraph.json]",
   "polishes a multiple sequence alignment graph",
   """zero or one pangraph file (native json)
      if no file, reads from stdin
      stream can be optionally gzipped.""",
   Arg[],
   (args) -> let
       path = parse(Polish, args)
       length(path) > 1

       graph = if path === nothing
           Graphs.unmarshal(stdin)
       elseif length(path) == 1
           path = path[1]
           !isfile(path) && error("file '$(path)' not found")
           open(Graphs.unmarshal, path)
       else
           usage(Polish)
           return 2
       end

       if !Graphs.havecommand("mafft")
           panic("external command mafft not found. please install before running polish step\n")
       end
       Graphs.realign!(graph)

       marshal(stdout, graph, :json)
       return 0
   end
)
