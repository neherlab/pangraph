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

def main(args):
    for arg in args:
        in_dir = f"data/{arg}/assemblies"
        if not os.path.exists(in_dir):
            print(f"{in_dir} doesn't exist. skipping...")
            continue

        out_dir = f"data/{arg}-plasmid/assemblies"
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)

        for path in glob(f"{in_dir}/*.f?a*"):
            with open(path, 'rt') as fd, open(f"{out_dir}/{os.path.basename(path).replace('.gz', '')}", 'w') as wtr:
                for i, rec in enumerate(parse_fasta(fd)):
                    if i == 0:
                        continue
                    wtr.write(str(rec))
                    wtr.write('\n')

if __name__ == "__main__":
    main(sys.argv[1:])
