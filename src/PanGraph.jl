module PanGraph

ENV["PYTHON"] = ""

using GZip
using Rematch
using Random: seed!

include("graph.jl")
using .Graphs

import .Graphs.Utility: read_fasta, write_fasta

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
        return main(ARGS)
    catch
        # TODO: more sophisticated error handling
        return 1
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(ARGS)
end

end
