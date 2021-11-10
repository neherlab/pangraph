Build = Command(
   "build",
   "pangraph build <options> [arguments]",
   "builds a multiple sequence alignment graph",
   """one or more fasta files.
      files can be optionally gzipped.
      multiple records within one file are treated as seperate genomes""",
   [
    Arg(
        Int,
        "minimum length",
        (short="-l", long="--len"),
        "minimum block size for alignment graph (in nucleotides)",
        100,
    ),
    Arg(
        Float64,
        "block junction cost",
        (short="-m", long="--mu"),
        "energy cost for introducing junction due to alignment merger",
        100,
    ),
    Arg(
        Float64,
        "block diversity cost",
        (short="-b", long="--beta"),
        "energy cost for interblock diversity due to alignment merger",
        20,
    ),
    Arg(
        Bool,
        "circular genomes",
        (short="-c", long="--circular"),
        "toggle if input genomes are circular",
        false,
    ),
    Arg(
        String,
        "distance calculator",
        (short="-d", long="--distance-backend"),
        """
        backend to use to estimate pairwise distance for guide tree
                    recognized options:
                        native
                        mash""",
        "native",
    ),
   ],

   (args) -> let
       files = parse(Build, args)
       files === nothing && return 2

       minblock = arg(Build, "-l")
       circular = arg(Build, "-c")
       μ = arg(Build, "-m")
       β = arg(Build, "-b")

       energy = (aln) -> let
           len = aln.length
           len < minblock && return Inf

           cuts(hit) = (hit.start > minblock) + ((hit.length-hit.stop) > minblock)

           ncuts = cuts(aln.qry)+cuts(aln.ref)
           nmuts = aln.divergence*aln.length

           return -len + μ*ncuts + β*nmuts
       end

       graph(io) = graphs(io; circular=circular)
       isolates  = (G for file in files for G ∈ (endswith(file,".gz") ? GZip.open(graph,file) : open(graph,file)))

       compare = @match arg(Build, "-d") begin
           "native" => Graphs.Mash.distance
           "mash"   => begin
               if !Graphs.havecommand("mash")
                   panic("external command mash not found. either install or use native backend\n")
               end
               Graphs.mash
           end
           _        => begin
                # XXX: hacky...
                Build.arg[5].value = "native" 
                usage(Build)
                exit(1)
           end
       end

       graph = Graphs.align(isolates...; compare=compare, energy=energy, minblock=minblock)
       marshal(stdout, graph; fmt=:json)
   end
)
