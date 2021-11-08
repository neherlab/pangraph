Polish = Command(
   "polish",
   "pangraph polish <options> [pangraph.json]",
   "polishes a multiple sequence alignment graph",
   """zero or one pangraph file (native json)
      if no file, reads from stdin
      stream can be optionally gzipped.""",
   Arg[],
   function(args)
       path = parse(Polish, args)
       length(path) > 1 && return 2

       graph = load(path, Polish)

       if !Shell.havecommand("mafft")
           panic("external command mafft not found. please install before running polish step\n")
       end
       Graphs.realign!(graph)

       marshal(stdout, graph, :json)
       return 0
   end
)
