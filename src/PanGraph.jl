module PanGraph

using Rematch

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

# ------------------------------------------------------------------------
# subcommands and arguments

include("args.jl")
using .Commands

# ---------------------------
# All subcommands go here

include("build.jl")
include("generate.jl")

pangraph = Command(
    "pangraph",
    "pangraph <command> [arguments]",
    "pangraph is a tool for aligning large sets of genomes in the presence of horizontal gene transfer",
    "passed directly to the chosen command",
    [
     Build,
     Generate,
    ],
)

function main(args)
    if length(args) == 0
        usage(pangraph)
        return 2
    end

    return run(pangraph, parse(pangraph, args))
end

function julia_main()::Cint
    try
        main(ARGS)
        return 0
    catch
        # TODO: more sophisticated error handling
        return 1
    end
end

end

if !isdefined(Base, :active_repl)
    PanGraph.main(ARGS)
end
