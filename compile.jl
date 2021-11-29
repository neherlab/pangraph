import Pkg; Pkg.add("PackageCompiler")
import PackageCompiler

PackageCompiler.create_app(
    ".",
    "pangraph",
    app_name="pangraph",
    force=true,
    filter_stdlibs=true,
    precompile_execution_file="trace.jl"
)
