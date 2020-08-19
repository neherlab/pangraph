#!/usr/bin/env python3
"""
script to filter plasmids and chromosomes from full genome assemblies
"""

import os
import sys
import gzip
import builtins

from glob import glob

sys.path.insert(0, os.path.abspath('.')) # gross hack
from pangraph.utils import parse_fasta, breakpoint

from Bio import SeqIO

def open(path, *args, **kwargs):
    if path.endswith('.gz'):
        return gzip.open(path, *args, **kwargs)
    else:
        return builtins.open(path, *args, **kwargs)

# def main(args):
#     for d in args:

from time import time
if __name__ == "__main__":
    t0 = time()
    for path in glob("data/staph/assemblies/*.fna.gz"):
        with open(path, 'rt') as fd:
            seqs = [record for record in SeqIO.parse(fd, 'fasta')]
    t1 = time()
    print(f"bio parser took {t1 - t0} seconds")

    t0 = time()
    for path in glob("data/staph/assemblies/*.fna.gz"):
        with open(path, 'rt') as fd:
            seqs = [record for record in parse_fasta(fd)]
    t1 = time()
    print(f"my parser took {t1 - t0} seconds")

    # main(sys.argv[1:])
