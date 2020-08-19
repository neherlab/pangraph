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

def open(path, *args, **kwargs):
    if path.endswith('.gz'):
        return gzip.open(path, *args, **kwargs)
    else:
        return builtins.open(path, *args, **kwargs)

# def main(args):
#     for d in args:

from time import time
if __name__ == "__main__":
    for path in glob("data/staph/assemblies/*.fna.gz"):
        with open(path, 'rt') as fd, open("test.fa", 'w') as wtr:
            for rec in parse_fasta(fd):
                # print(str(rec))
                wtr.write(str(rec))
                wtr.write('\n')
            seqs = [record for record in parse_fasta(fd)]
        break

    # main(sys.argv[1:])
