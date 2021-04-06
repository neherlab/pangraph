Generate = Command(
   "generate",
   "pangraph generate <options>",
   "generates a simulated sequence alignment graph from a neutral WF model with HGT",
   "",
   [
    Arg(
        Int,
        "population size",
        (short="-n", long="--size"),
        "simulated population size",
        100,
    ),
    Arg(
        Float64,
        "genome length",
        (short="-l", long="--len"),
        "average (and initial) simulated genome length",
        100,
    ),
    Arg(
        Float64,
        "mutation rate",
        (short="-m", long="--mu"),
        "rate of mutations per site per genome per generation",
        1e-5,
    ),
    Arg(
        Int,
        "generations simulation",
        (short="-t", long="--time"),
        "number of generations simulated under WF model",
        35,
    ),
   ]
)
