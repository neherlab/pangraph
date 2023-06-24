Build = Command(
    "build",
    "pangraph build <options> [arguments]",
    "align genomes into a multiple sequence alignment graph",
    """zero, one or more fasta files.
       files can be optionally gzipped.
       multiple records within one file are treated as seperate genomes.
       if no file is given, reads fasta record from stdin""",
    [
        Arg(
            Int,
            "minimum length",
            (short = "-l", long = "--len"),
            "minimum block size for alignment graph (in nucleotides)",
            100,
        ),
        Arg(
            Float64,
            "block junction cost",
            (short = "-a", long = "--alpha"),
            "energy cost for introducing junction due to alignment merger",
            100,
        ),
        Arg(
            Float64,
            "block diversity cost",
            (short = "-b", long = "--beta"),
            "energy cost for interblock diversity due to alignment merger",
            10,
        ),
        Arg(
            Bool,
            "circular genomes",
            (short = "-c", long = "--circular"),
            "toggle if input genomes are circular",
            false,
        ),
        Arg(
            Bool,
            "enforce uppercase",
            (short = "-u", long = "--upper-case"),
            "transforms all sequence to upper case",
            false,
        ),
        Arg(
            String,
            "pairwise sensitivity",
            (short = "-s", long = "--sensitivity"),
            "used to set pairwise alignment sensitivity\n\trecognized options: [5, 10, 20]",
            "10",
        ),
        Arg(
            Int,
            "maximum self maps",
            (short = "-x", long = "--max-self-map"),
            "maximum number of self mappings to consider per pairwise graph merger",
            100,
        ),
        Arg(
            String,
            "distance calculator",
            (short = "-d", long = "--distance-backend"),
            "backend to use to estimate pairwise distance for guide tree\n\trecognized options: [native, mash]",
            "native",
        ),
        Arg(
            String,
            "alignment kernel",
            (short = "-k", long = "--alignment-kernel"),
            "backend to use for pairwise genome alignment\n\trecognized options: [minimap2, mmseqs]",
            "minimap2",
        ),
        Arg(
            Int,
            "k-mer length",
            (short = "-K", long = "--kmer-length"),
            "kmer length, only used for mmseqs2 alignment kernel. If not specified will use mmseqs default.",
            0,
        ),
        Arg(
            Bool,
            "verify that input genomes can be reconstructed from output",
            (short = "-t", long = "--test"),
            "toggle to activate consistency check at each step of the graph merging process",
            false,
        ),
    ],
    (args) -> let
        files = parse(Build, args)
        files = if files === nothing || length(files) == 0
            ["/dev/stdin"]
        else
            files
        end

        minblock = arg(Build, "-l")
        circular = arg(Build, "-c")
        uppercase = arg(Build, "-u")

        α = arg(Build, "-a")
        β = arg(Build, "-b")

        energy = function (aln)
            len = aln.length
            len < minblock && return Inf

            cuts(hit) = (hit.start > minblock) + ((hit.length - hit.stop) > minblock)

            ncuts = cuts(aln.qry) + cuts(aln.ref)
            nmuts = aln.divergence * aln.length

            return -len + α * ncuts + β * nmuts
        end

        maxiter = arg(Build, "-x")
        sensitivity = @match arg(Build, "-s") begin
            "5" => "asm5"
            "10" => "asm10"
            "20" => "asm20"
            _ => begin
                usage(Build)
                exit(1)
            end
        end

        graph(io) = graphs(io; circular = circular, upper = uppercase)
        isolates = (G for file in files for G ∈ open(graph, file))

        compare = @match arg(Build, "-d") begin
            "native" => Graphs.Mash.distance
            "mash" => begin
                if !Graphs.havecommand("mash")
                    panic(
                        "external command mash not found. either install or use native backend\n",
                    )
                end
                Graphs.mash
            end
            _ => begin
                # XXX: hacky...
                Build.arg[8].value = "native"

                usage(Build)
                exit(1)
            end
        end

        aligner(contigs₁, contigs₂) = @match arg(Build, "-k") begin
            "minimap2" => Minimap.align(contigs₁, contigs₂, minblock, sensitivity)
            "mmseqs" => let
                if !Shell.havecommand("mmseqs")
                    panic(
                        "external command mmseqs not found. please install before running build step with mmseqs backend\n",
                    )
                end
                kmer_len = arg(Build, "-K")
                MMseqs.align(contigs₁, contigs₂, kmer_len)
            end
            _ => error("unrecognized alignment kernel")
        end

        # if testing is set to true, collect dictionary of reference sequences for comparisons
        reference =
            arg(Build, "-t") ?
            Dict(begin
                    key = collect(keys(G.sequence))[1]
                    seq = sequence(collect(values(G.sequence))[1])
                    key => seq
                end for G in isolates) : nothing

        graph = Graphs.align(
            aligner,
            isolates...;
            compare = compare,
            energy = energy,
            minblock = minblock,
            maxiter = maxiter,
            reference = reference,
            verbose = false,
        )
        finalize!(graph)

        # test graph consistency
        arg(Build, "-t") && Graphs.consistency_check(graph)

        marshal(stdout, graph; fmt = :json)
    end,
)
