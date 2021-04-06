import Pkg
Pkg.add("PackageCompiler")
import PackageCompiler

PackageCompiler.create_app(
    ".",
    ".build", 
    app_name="pangraph",
    filter_stdlibs=true, 
    precompile_execution_file="precompile.jl"
)
