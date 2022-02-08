import Pkg; Pkg.add("PackageCompiler")
import PackageCompiler

PackageCompiler.create_app(
    ".",
    "pangraph",
    executables=["pangraph"=>"julia_main"],
    force=true,
    include_lazy_artifacts=false,
    precompile_execution_file="trace.jl"
)
