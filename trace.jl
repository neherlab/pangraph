import PanGraph

home = ENV["HOME"]
root = "$home/src/pangraph/data/synthetic"

# usage
PanGraph.main([])              # main usage
PanGraph.main(["build"])       # build usage
PanGraph.main(["export"])      # export usage
PanGraph.main(["generate"])    # generate usage
PanGraph.main(["marginalize"]) # marginalize usage
PanGraph.main(["polish"])      # polish usage

# build
PanGraph.main(["build", "--circular", "-b", "0", "-m", "0", "$root/test.fa"])

# export
PanGraph.main(["export", "-o", "$root/export", "$root/test.json"])

# export
PanGraph.main(["generate", "-o", "$root/export.json", "$root/test.fa"])

# marginalize
PanGraph.main(["marginalize", "-o", "$root/pairs", "$root/test.json"])

# polish
PanGraph.main(["polish", "$root/test.json"])
