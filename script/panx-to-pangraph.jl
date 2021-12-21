module PanX

using Printf
using JSON, GZip
using PanGraph

struct Gene
    id     :: Int
    tag    :: String
    low    :: Int
    high   :: Int
    strand :: Bool
end

function find(database, tag)
    for (i,cluster) in enumerate(database)
        tag ∈ cluster && return i
    end
    return nothing
end

askey(name,tag) = "$(name)|$(tag)"

function gene(line, name, database)
    elt = split(line)

    strand = elt[1] == "+"
    low, high = parse.(Int, split(elt[2],".."))
    tag = elt[3]
    id  = find(database,askey(name,tag))
    id !== nothing || return nothing

    return Gene(id, tag, low, high, strand)
end

struct Genome
    name     :: String
    sequence :: Array{UInt8,1}
    gene     :: Array{Gene,1}
end

Genome(name::String, sequence::Array{UInt8,1}) = Genome(name, sequence, Gene[])

function database(path)
    data = Set{String}[]
    meta = String[]
    re   = r"(N[A-Z]_[A-Z0-9]+)_"
    sub  = s"\1|"
    open(path) do io
        genes = JSON.parse(io)
        for (i,gene) in enumerate(genes)
            loci = map(split(gene["locus"])) do g
                replace(g, re => sub)
            end
            push!(data, Set(loci))
            push!(meta, gene["msa"])
        end
    end
    return (
        loci=data,
        meta=meta,
        path=dirname(path),
    )
end

function alignment(io::IO)
    record = [ (name=rec.name, seq=rec.seq) for rec in PanGraph.Utility.read_fasta(io)]
    return hcat(map(r->r.seq,record)...), map(record) do r
        entry = split(r.name, "-")
        return "$(entry[1])|$(entry[2])"
    end
end

function pangraph(genomes, genes)
    vertex = [
        let
            path = "$(genes.path)/geneCluster/$(base)_na_aln.fa.gz"
            aln, isos = GZip.open(alignment, path)

            ref = PanGraph.Utility.make_consensus(aln)
            blk = PanGraph.Graphs.Block(ref)
            iso = [ PanGraph.Graphs.Node(blk, true) for _ in isos ]

            blk.uuid = @sprintf("Gene_%05d", i)
            blk.gaps,
            blk.mutate,
            blk.delete,
            blk.insert,
            blk.sequence = PanGraph.Utility.alignment_alleles(ref, aln, iso)

            (block=blk, node=Dict( name => node for (name,node) in zip(isos,iso)))
        end for (i,(gene,base)) in enumerate(zip(genes.loci,genes.meta))
    ]

    path = [
        let
            nodes = [
             let
                 node = pop!(vertex[gene.id].node, askey(genome.name, gene.tag))
                 node.strand = gene.strand
                 node
             end for gene in genome.gene
            ]
            breaks = [ x for gene in genome.gene for x in [gene.low, gene.high] ]
            PanGraph.Graphs.Path(genome.name, nodes, 0, true, breaks)
        end for genome in values(genomes)
    ]

    # remove any unpaired nodes we initially added
    # delete empty blocks
    index = Int[]
    for (i,v) in enumerate(vertex)
        for n in values(v.node)
            pop!(v.block, n)
        end

        if length(v.block.mutate) == 0
            push!(index, i)
        end
    end
    deleteat!(vertex, index)

    block = map(v->v.block, vertex)
    graph = PanGraph.Graph(Dict(b.uuid => b for b in block), Dict(p.name => p for p in path) )

    return graph
end

function main(database)
    genomes = Dict{String,Genome}()
    while !eof(stdin)
        line = readline(stdin)
        line[1] == '>' || error("invalid stream detected")

        name = line[2:end]
        sbuf = IOBuffer()

        # grab sequence
        @label moreseq
        line = readline(stdin)
        if line != "ENDSEQ"
            print(sbuf, line)
            @goto moreseq
        end
        genome = Genome(name, take!(sbuf))

        # grab loci
        @label moregene
        line = readline(stdin)
        if line != "ENDGENE"
            g = gene(line, name, database.loci)
            g !== nothing && push!(genome.gene, g)
            @goto moregene
        end

        genomes[name] = genome
    end

    return pangraph(genomes, database)
end

# NOTE: just for debugging
function echo()
    while !eof(stdin)
        println(readline(stdin))
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    database(ARGS[1]) |> main |> (G) -> PanGraph.Graphs.marshal(stdout, G; fmt=:json)
end

end
