import os, sys
import argparse
import pickle
import pyfaidx as fai
import numpy as np

from glob  import glob

from py.graph import Graph
from py.kmers import Tree, parse

def empirical():
    seqs = fai.Fasta("data/all_plasmids_filtered.fa")
    # M, nms = parse("data/kmerdist.txt")
    # T      = Tree.nj(M, nms)
    # T.align(seqs)
    for f in glob("data/graph/nwk/*.nwk"):
        if "024" not in f:
            continue

        betas = [1, 2, 5, 10, 25, 50, 75, 100, 250, 500, 750, 1000]
        mus   = [1e0, 5e0, 1e1, 5e1, 1e2, 5e2, 1e3, 5e3, 1e4]

        for mu in mus:
            for beta in betas:
                name = os.path.basename(f).split(".")[0]
                print(f"--> Building graph for mu={mu}, beta={beta}")
                g, ok = Graph.fromnwk(f, seqs, save=False, verbose=False, mu=mu, beta=beta)

                g.name = f"{name}_m{int(mu):d}_b{int(beta):d}"
                g.tojson()
                g.tofasta()
                pickle.dump(g.todict(), open(f"data/graph/pkl/{g.name}.pkl", "wb"))

        # if ok:
            # g.compilesuffix()
            # g.computepairdists()

        # print(f"-->Saving graph to pickle")
        # print(f"-->Saving graph to json")

def synthetic():
    f    = "test.nwk"
    seqs = fai.Fasta("test.fasta")

    g, ok = Graph.fromnwk(f, seqs)
    g.name = "test"
    g.tofasta()
    pickle.dump(g.todict(), open(f"{g.name}.pkl", "wb"))

    return g

if __name__ == "__main__":
    g = synthetic()
