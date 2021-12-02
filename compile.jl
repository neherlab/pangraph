import Pkg; Pkg.add("PackageCompiler")
import PackageCompiler

PackageCompiler.create_app(
    ".",
    "pangraph",
    app_name="pangraph",
    force=true,
    precompile_execution_file="trace.jl"
)
