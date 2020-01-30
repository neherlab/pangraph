import os, sys
import argparse
import pickle
import pyfaidx as fai
import numpy as np

from glob  import glob

from py.graph import Graph
from py.kmers import Tree, parse

if __name__ == "__main__":
    seqs = fai.Fasta("data/all_plasmids_filtered.fa")
    # M, nms = parse("data/kmerdist.txt")
    # T      = Tree.nj(M, nms)
    # T.align(seqs)
    for f in glob("data/graph/nwk/*.nwk"):
        name = os.path.basename(f).split(".")[0]
        print(f"Processing {name}")
        print(f"-->Building graph")
        g, ok = Graph.fromnwk(f, seqs, save=False, verbose=False)
        if ok:
            g.compilesuffix()
            g.computepairdists()

        print(f"-->Saving graph to pickle")
        pickle.dump(g.todict(), open(f"data/graph/pkl/{name}.pkl", "wb"))
        print(f"-->Saving graph to json")
        g.name = name
        g.tojson()
        print(f"-->Saving graph to fasta")
        g.tofasta()

