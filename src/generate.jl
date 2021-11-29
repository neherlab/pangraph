Generate = Command(
   "generate",
   "pangraph generate <options> [input fasta file]",
   "generates a simulated sequence alignment graph",
   """zero or one fasta file.
      if zero, sequences are assumed to be supplied via stdin.""",
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
        "mutation rate",
        (short="-m", long="--snp-rate"),
        "rate of mutations per site per genome per generation",
        1e-5,
    ),
    Arg(
        Float64,
        "HGT rate",
        (short="-r", long="--hgt-rate"),
        "rate of horizontal transfer events per genome per generation",
        0,
    ),
    Arg(
        Float64,
        "Deletion rate",
        (short="-d", long="--delete-rate"),
        "rate of deletion events per genome per generation",
        0,
    ),
    Arg(
        Float64,
        "Inversion rate",
        (short="-i", long="--invert-rate"),
        "rate of inversion events per genome per generation",
        0,
    ),
    Arg(
        String,
        "Graph output path",
        (short="-o", long="--graph-output"),
        "path where to output graph",
        "",
    ),
    Arg(
        Int,
        "Generations to simulate",
        (short="-t", long="--time"),
        "number of generations simulated under WF model",
        35,
    ),
   ],
   function(args)
        if args === nothing || length(args) == 0 
            usage(Generate)
            return 2
        end

        input = parse(Generate, args)
        input = if (input === nothing) 
            stdin
        elseif length(input) == 1
            input[1]
        else
            return 2
        end

        ancestors = open(input) do io
            map(r->r.seq, read_fasta(io))
        end

        N = length(ancestors)
        L = Int(round(sum(length.(ancestors)) / N))
        σ = L / 10

        evolve! = Simulation.model(
            Simulation.Params(;
                N   = N,
                L   = L,
                σₗ  = σ,
                snp = arg(Generate,"-m"),
                hgt = arg(Generate,"-r"),
                del = arg(Generate,"-d"),
                inv = arg(Generate,"-i"),
            )
        )

        T    = arg(Generate, "-t")
        path = arg(Generate, "-o")

        sequences, graph = Simulation.run(evolve!, T, ancestors; graph=path!="")

        for (i, sequence) in enumerate(sequences)
            write_fasta(stdout, "isolate_$(i)", sequence)
        end

        if length(path) > 0
            open(path,"w") do io
                marshal(io, graph; fmt=:json)
            end
        end
    end
)
