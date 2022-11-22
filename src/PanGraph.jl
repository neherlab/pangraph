module PanGraph

ENV["PYTHON"] = ""

using GZip
using Rematch
using Random: seed!

# ------------------------------------------------------------------------
# types

Maybe{T} = Union{Nothing,T}

"""
	struct PanContigs
		name     :: T
		sequence :: T
	end

A synonym for a consensus sequence of `Block`.
"""
struct PanContigs{T <: AbstractArray{S} where S <: AbstractString}
    name     :: T
    sequence :: T
end

# paf alignment pair
"""
	mutable struct Hit
		name::String
		length::Int
		start::Int
		stop::Int
		seq::Maybe{Array{UInt8,1}}
	end

Hit is one side of a pairwise alignment between homologous sequences.
"""
mutable struct Hit
    name::String
    length::Int
    start::Int
    stop::Int
    seq::Maybe{Array{UInt8,1}}
end

"""
	mutable struct Alignment{T <: Union{String,Nothing,Array{Tuple{Int,Char}}}}
		qry::Hit
		ref::Hit
		matches::Int
		length::Int
		quality::Int
		orientation::Bool
		cigar::T
		divergence::Union{Float64,Nothing}
		align::Union{Float64,Nothing}
	end

Alignment is a pairwise homologous alignment between two sequences.
"""
mutable struct Alignment
    qry::Hit
    ref::Hit
    matches::Int
    length::Int
    quality::Int
    orientation::Bool
    cigar::Union{String,Nothing,Array{Tuple{Int,Char}}}
    divergence::Union{Float64,Nothing}
    align::Union{Float64,Nothing}
end


# ------------------------------------------------------------------------
# local imports

include("graph.jl")
using .Graphs

import .Graphs.Utility: read_fasta, write_fasta

# alignment kernels
include("minimap.jl")
include("mmseqs.jl")

include("simulate.jl")
using .Simulation

export Graphs, Simulation

# ------------------------------------------------------------------------
# errors

struct Error <: Exception
    msg::AbstractString
end
Base.showerror(io::IO, e::Error) = print(io, "PanGraph Error: ", e.msg)

function panic(msg...)
    print(stderr, string(msg...))
    exit(2)
end

function open(func, path, args...)
    endswith(path, ".gz") && return GZip.open(func, path, args...)
    return Base.open(func, path, args...)
end

function load(path, cmd)
   if path === nothing
       unmarshal(stdin)
   elseif length(path) == 1
       path = path[1]
       !isfile(path) && error("file '$(path)' not found")
       open(unmarshal, path)
   else
       usage(cmd)
       return 2
   end
end

# ------------------------------------------------------------------------
# subcommands and arguments

include("args.jl")
using .Commands

# ---------------------------
# All subcommands go here

include("build.jl")
include("generate.jl")
include("help.jl")
include("version.jl")
include("polish.jl")
include("marginalize.jl")
include("export.jl")

Dispatch = Command(
    "pangraph",
    "pangraph <command> [arguments]",
    "a tool for aligning large sets of closely related genomes in the presence of horizontal gene transfer",
    "passed directly to the chosen command",
    [
     Build,
     Generate,
	 Help,
     Version,
     Polish,
     Marginalize,
     Export,
    ],
)

function main(args)
    if length(args) == 0
        usage(Dispatch)
        return 2
    end

    return run(Dispatch, parse(Dispatch, args))
end

function julia_main()::Cint
    try
        main(ARGS)
    catch
        Base.invokelatest(Base.display_error, Base.catch_stack())
        return 1
    end

    return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(ARGS)
end

end
