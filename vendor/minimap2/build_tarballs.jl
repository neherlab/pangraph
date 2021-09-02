using BinaryBuilder

name = "minimap2"
version = v"2.22"

sources = [
    GitSource("https://github.com/lh3/minimap2", "7358a1ead1adfa89a2d3d0e72ffddd05732f9850"),
    DirectorySource("./static")
]

script = raw"""
cd $WORKSPACE/srcdir/minimap2
sed 's/SO/${dlext}/g' $WORKSPACE/srcdir/patches/Makefile.diff | patch Makefile 
INC=${includedir} make libminimap2.${dlext}
cp libminimap2.${dlext} ${prefix}/lib
"""

platforms = filter(supported_platforms()) do platform
    (Sys.islinux(platform) | Sys.isapple(platform)) && platform.tags["arch"] == "x86_64"
end

products = [
    LibraryProduct("libminimap2", :libminimap2),
]

dependencies = [
    Dependency("Zlib_jll"),
]

build_tarballs(ARGS, name, version, sources, script, platforms, products, dependencies)
