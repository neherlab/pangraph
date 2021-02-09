Generate = Command(
   "generate",
   "pangraph generate <options>",
   "generates a simulated sequence alignment graph from a neutral WF model with HGT",
   "",
   [
    Arg(
        Int,
        "population size",
        (long="-n", short="--size"),
        "simulated population size",
        100,
    ),
    Arg(
        Float64,
        "genome length",
        (long="-l", short="--len"),
        "average (and initial) simulated genome length",
        100,
    ),
    Arg(
        Float64,
        "mutation rate",
        (long="-m", short="--mu"),
        "rate of mutations per site per genome per generation",
        1e-5,
    ),
    Arg(
        Int,
        "generations simulation",
        (long="-t", short="--time"),
        "number of generations simulated under WF model",
        35,
    ),
   ]
)
