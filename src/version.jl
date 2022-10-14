# using Pkg
# const PANGRAPH_VERSION = begin
#     dir = dirname(string(first(methods(PanGraph.eval)).file))
#     project_file = joinpath(dir, "..", "Project.toml")
#     Pkg.TOML.parsefile(project_file)["version"]
# end

# using Pkg
# pangraph_pkg_status = sprint(io -> Pkg.status("PanGraph"; io = io))
# re = r"PanGraph (\S+)"
# const PANGRAPH_VERSION = match(re, pangraph_pkg_status)[1]

const PANGRAPH_VERSION = "v0.6.0"

Version = Command(
    "version",
    "pangraph version",
    "prints pangraph version on stderr",
    "none",
    Arg[],
    function (args)
        command = parse(Version, args)
        command !== nothing && return usage(Version)
        println(stderr, "Pangraph version $PANGRAPH_VERSION")
    end,
)

