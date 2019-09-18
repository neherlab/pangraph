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
# Global constants/variables

qryfa = "tmp/qry.fasta"

# ------------------------------------------------------------------------
# Functions

def aln(qry, ref, out):
    os.system(f"minimap2 -x asm5 -D -c {qry} {ref} 1>{out}.paf 2>log")

def writefa(seq, out):
    with open(out, 'w') as fh:
        fh.write(f">{str(seq.name)}\n")
        fh.write(str(seq[:]))

def main(args):
    # Fasta file must be indexed for random access
    seqs = fai.Fasta(args.fa)
    for isolate in seqs.keys():
        if os.path.isfile(f"{args.dir}/{isolate}.paf"):
            continue
        print(isolate)
        writefa(seqs[isolate], qryfa)
        aln(qryfa, args.fa, f"{args.dir}/{isolate}")


# ------------------------------------------------------------------------
# Main point of entry

parser = argparse.ArgumentParser(description = "", usage = "map each isolate against entire database")
parser.add_argument("fa",  type = str, help = "fasta file")
parser.add_argument("dir", type = str, help = "output directory")

if __name__ == "__main__":
    args = parser.parse_args()
    main(args)
