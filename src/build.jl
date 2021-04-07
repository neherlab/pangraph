using GZip

include("graph.jl")
using .Graphs

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
   ],

   (args) -> begin
       files  = parse(Build, args)

       maxgap = Build.arg[1].value
       energy = (aln) -> let
           len = aln.length
           len < maxgap && return Inf

           cuts(hit) = (hit.start > maxgap) + ((hit.length-hit.stop) > maxgap)

           ncuts = cuts(aln.qry)+cuts(aln.ref)
           nmuts = aln.divergence*aln.length

           return -len + Build.arg[2].value*ncuts + Build.arg[3].value*nmuts
       end
       isolates = (G for file in files for G ∈ (endswith(file,".gz") ? GZip.open(graphs,file) : open(graphs,file)))

       graph = Graphs.align(isolates...; energy=energy, maxgap=maxgap)
       # marshal(stdout, graph, :json)
   end
)
