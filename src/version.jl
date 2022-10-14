using Pkg

# extract pangraph version using pkg status
const PANGRAPH_VERSION = begin
    pangraph_pkg_status = sprint(io -> Pkg.status("PanGraph"; io = io))
    re = r"PanGraph (\S+)"
    String(match(re, pangraph_pkg_status)[1])
end

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