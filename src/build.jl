include("graph.jl")
import .Graphs

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
        (long="-l", short="--len"),
        "minimum block size for alignment graph (in nucleotides)",
        50,
    ),
    Arg(
        Float64,
        "block junction cost",
        (long="-m", short="--mu"),
        "energy cost for introducing junction due to alignment merger",
        100,
    ),
    Arg(
        Float64,
        "block diversity cost",
        (long="-b", short="--beta"),
        "energy cost for interblock diversity due to alignment merger",
        20,
    ),
    Arg(
        Bool,
        "circular genomes",
        (long="-c", short="--circular"),
        "toggle if input genomes are circular",
        false,
    ),
   ],

   (args) -> begin
       files = parse(Build, args)
       Graphs.test()
   end
)