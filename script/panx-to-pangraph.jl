module PanX

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
    sequence :: Array{UInt8,1}
    gene     :: Array{Gene,1}
end
Genome(sequence::Array{UInt8,1}) = Genome(sequence, Gene[])

function main(path)
    database = Set{String}[]
    open(path) do io
        for (i,gene) in enumerate(eachline(path))
            push!(database, Set(split(gene)))
        end
    end
    @show database

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
        genome = Genome(take!(sbuf))

        # grab loci
        @label moregene
        line = readline(stdin)
        if line != "ENDGENE"
            g = gene(line, name, database)
            g !== nothing && push!(genome.gene, g)
            @goto moregene
        end

        genomes[name] = genome
    end

    @show genomes

    # TODO: turn into a pangraph
end

# NOTE: just for debugging
function echo()
    while !eof(stdin)
        println(readline(stdin))
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(ARGS[1])
    # echo()
end

end
