import PanGraph

cdir = pwd()
home = replace(cdir,"/pangraph/lib/julia"=>"")
root = "$home/data/synthetic"

# usage
PanGraph.main([])                      # main usage
PanGraph.main(["help", "build"])       # build usage
PanGraph.main(["help", "export"])      # export usage
PanGraph.main(["help", "generate"])    # generate usage
PanGraph.main(["help", "marginalize"]) # marginalize usage
PanGraph.main(["help", "polish"])      # polish usage
PanGraph.main(["help", "version"])		   # version usage
PanGraph.main(["help", "help"])		   # help usage

# build (native - mmseqs)
PanGraph.main(["build", "-v", "-c", "-u", "-b", "0", "-a", "0", "$root/test.fa"])
# PanGraph.main(["build", "-c", "-u", "-d", "mash", "$root/test.fa"])
PanGraph.main(["build", "-v", "-c", "-u", "-b", "0", "-a", "0", "-k", "mmseqs", "-K", "8", "$root/test.fa"])

# export
PanGraph.main(["export", "-o", "$root/export", "$root/test.json"])

# generate
PanGraph.main(["generate", "-o", "$root/export.json", "$root/test.fa"])

# marginalize
PanGraph.main(["marginalize", "-v", "-o", "$root/pairs", "$root/test.json"])

# polish
PanGraph.main(["polish", "-c", "-l", "10000", "$root/test.json"])