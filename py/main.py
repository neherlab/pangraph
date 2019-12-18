import argparse
import pyfaidx as fai

from glob  import glob
from graph import Graph

# TODO: Parse in user given directory.

nwkdir = f"data/graph/nwk"
if __name__ == "__main__":
    seqs = fai.Fasta("data/seq/all.fasta")
    for path in glob(f"{nwkdir}/*.nwk"):
        print(f"Analyzing {path}")
        g = Graph.fromnwk(path, seqs, save=True, verbose=False)
        # Graph.cleanbld()
