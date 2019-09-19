import os
import argparse, gzip
from enum import IntEnum
from collections import defaultdict, OrderedDict

import numpy as np
import matplotlib.pylab as plt
import seaborn as sns

import pyfaidx as fai

from graph import Graph
from cigar import Cigar
from util  import parse_paf
from interval import Partition

# ------------------------------------------------------------------------
# Functions

def writefa(seq, out):
    with open(out, 'w') as fh:
        fh.write(f">{str(seq.name)}\n")
        fh.write(str(seq[:]))

def main(args):
    # Fasta file must be indexed for random access
    seqs = fai.Fasta(args.fa)
    for isolate in seqs.keys():
        outpath = f"{args.dir}/{isolate}.fasta"
        if os.path.isfile(outpath):
            continue
        writefa(seqs[isolate], outpath)

# ------------------------------------------------------------------------
# Main point of entry

parser = argparse.ArgumentParser(description = "", usage = "map each isolate against entire database")
parser.add_argument("fa",  type = str, help = "fasta file")
parser.add_argument("dir", type = str, help = "output directory")

if __name__ == "__main__":
    args = parser.parse_args()
    main(args)
