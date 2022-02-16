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
PanGraph.main(["help", "help"])		   # help usage

# build
PanGraph.main(["build", "-c", "-u", "-b", "0", "-a", "0", "$root/test.fa"])

# export
PanGraph.main(["export", "-o", "$root/export", "$root/test.json"])

# export
PanGraph.main(["generate", "-o", "$root/export.json", "$root/test.fa"])

# marginalize
PanGraph.main(["marginalize", "-o", "$root/pairs", "$root/test.json"])

# polish
PanGraph.main(["polish", "-c", "-l", "10000", "$root/test.json"])

# help
PanGraph.main(["help", "build"])
