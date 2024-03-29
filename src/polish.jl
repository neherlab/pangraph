Polish = Command(
   "polish",
   "pangraph polish <options> [pangraph.json]",
   "realigns pancontigs of multiple sequence alignment graph",
   """zero or one pangraph file (native json)
      if no file, reads from stdin
      stream can be optionally gzipped.""",
   [
    Arg(
        Int,
        "maximum length",
        (short="-l", long="--length"),
        "cutoff above which we won't realign",
        typemax(Int),
    ),
    Arg(
        Bool,
        "preserve case",
        (short="-c", long="--preserve-case"),
        "ensure case (upper/lower) is preserved after realignment",
        false,
    ),
   ],
   function(args)
       path = parse(Polish, args)
       path = if (path === nothing || length(path) == 0)
           nothing
       elseif length(path) == 1
           path
       else
           return 2
       end

       graph = load(path, Polish)
       if !Shell.havecommand("mafft")
           panic("external command mafft not found. please install before running polish step\n")
       end

       case = arg(Polish, "-c")
       accept = function(blk)
           length(blk) ≤ arg(Polish, "-l") && Graphs.depth(blk) > 1
       end
       Graphs.realign!(graph; accept=accept, case=case)

       marshal(stdout, graph; fmt=:json)
       return 0
   end
)
