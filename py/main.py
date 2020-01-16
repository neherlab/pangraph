import argparse
import pyfaidx as fai

from glob  import glob
from graph import Graph
from kmers import Tree, parse

nwkdir = f"data/graph/nwk"
if __name__ == "__main__":
    seqs   = fai.Fasta("data/all_plasmids_filtered.fa")
    M, nms = parse("data/kmerdist.txt")
    T      = Tree.nj(M, nms)
    T.align(seqs)
    # g = Graph.fromnwk("data/kmer.nwk", seqs, save=False, verbose=False)
