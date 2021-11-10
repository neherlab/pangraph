include("panX.jl")

Export = Command(
   "export",
   "pangraph export <options> [arguments]",
   "exports a pangraph to a chosen file format(s)",
   """zero or one pangraph file (native json)
      if no file given, reads from stdin
      stream can be optionally gzipped.""",
   [
    Arg(
        Int,
        "minimum block length for edge connection",
        (short="-ell", long="--edge-minimum-length"),
        "blocks below this length cutoff will be ignored for edges in graph",
        200,
    ),
    Arg(
        Int,
        "maximum block length for edge connection",
        (short="-elu", long="--edge-maximum-length"),
        "blocks above this length cutoff will be ignored for edges in graph",
        typemax(Int),
    ),
    Arg(
        Int,
        "minimum block depth for edge connection",
        (short="-edl", long="--edge-mininum-depth"),
        "blocks below this depth cutoff will be ignored for edges in graph",
        0,
    ),
    Arg(
        Int,
        "maximum block depth for edge connection",
        (short="-edu", long="--edge-maximum-depth"),
        "blocks above this depth cutoff will be ignored for edges in graph",
        typemax(Int),
    ),
    Arg(
        Int,
        "minimum block length for export",
        (short="-ll", long="--minimum-length"),
        "blocks below this length cutoff will not be exported",
        0,
    ),
    Arg(
        Int,
        "maximum block length for export",
        (short="-lu", long="--maximum-length"),
        "blocks above this length cutoff will not be exported",
        typemax(Int),
    ),
    Arg(
        Int,
        "minimum block depth for export",
        (short="-dl", long="--mininum-depth"),
        "blocks below this depth cutoff will not be exported",
        0,
    ),
    Arg(
        Int,
        "maximum block depth for export",
        (short="-du", long="--maximum-depth"),
        "blocks above this depth cutoff will not be exported",
        typemax(Int),
    ),
    Arg(
        String,
        "output directory",
        (short="-o", long="--output-directory"),
        "relative path to directory where output will be stored",
        "export"
    ),
    Arg(
        String,
        "prefix",
        (short="-p", long="--prefix"),
        "basename of files",
        "pangraph"
    ),
    Arg(
        Bool,
        "no GFA export",
        (short="-ng", long="--no-export-gfa"),
        "do not emit GFA file",
        false,
    ),
    Arg(
        Bool,
        "export panX visualization",
        (short="-pX", long="--export-panX"),
        "emit vis directory to input to panX-visualization",
        false,
    ),
   ],

   function(args)
       path = parse(Export, args)
       path === nothing && return 2
       length(path) > 1 && return 2

       graph = load(path, Export)

       function cutoffs(prefix)
           (
              depth = (
                   min = arg(Export, "-$(prefix)dl"),
                   max = arg(Export, "-$(prefix)du"),
              ),
              length = (
                   min = arg(Export, "-$(prefix)ll"),
                   max = arg(Export, "-$(prefix)lu"),
              )
           )
       end

       node = cutoffs("")
       edge = cutoffs("e")

       filter = (
            connect = function(node)
                return (!(edge.depth.min  ≤ Blocks.depth(node.block)  ≤ edge.depth.max)
                     || !(edge.length.min ≤ Blocks.length(node.block) ≤ edge.length.max))
            end,
            output = function(segment)
                return (!(node.depth.min  ≤ segment.depth  ≤ node.depth.max)
                     || !(node.length.min ≤ length(segment.sequence) ≤ node.length.max))
            end
       )

       # Create directory if it doesn't exist
       directory = arg(Export, "-o")
       if !isdir(directory)
           mkpath(directory)
       end
       prefix = arg(Export, "-p")

       # GFA export (default)
       if !arg(Export, "-ng")
           Base.open("$(directory)/$(prefix).gfa", "w") do io
               marshal(io, graph; fmt=:gfa, opt=filter)
           end
       end

       # panX export (doesn't fit into marshal paradigm)
       if arg(Export, "-pX")
           if !Shell.havecommand("fasttree")
               panic("external command fasttree not found. please install before attempting to export as panX visualization\n")
           end

           vis = "$(directory)/vis"
           if !isdir(vis)
               mkpath(vis)
           end

           PanX.emit(graph, vis)
       end

        # TODO: others?
   end
)
