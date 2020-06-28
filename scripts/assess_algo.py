#!/usr/bin/env python3

import os
import sys
from io import StringIO

# gross hack
sys.path.insert(0, os.path.abspath('.'))

import subprocess

from glob import glob
from sys  import argv, exit
from Bio  import SeqIO
from pangraph.utils import parse_paf

argv0 = None

def usage():
    print(f"usage: {argv0} [directory]", file=sys.stderr)
    return 1

def main(args):
    if len(args) != 1 or not os.path.isdir(args[0]):
        exit(usage())
    dir      = args[0].rstrip("/")
    graphs   = glob(f"{dir}/graph_???.fa")
    anc_blks = f"{dir}/ancestral.fa"
    if len(graphs) < 1:
        exit(f"directory {dir}: missing aligned pangraphs")
    if not os.path.exists(anc_blks):
        exit(f"directory {dir}: missing ancestral block sequences")

    results = []
    for g in graphs:
        results.append(subprocess.run(
            ["minimap2", "-x", "asm5", str(g), str(anc_blks)],
            capture_output=True)
        )

    for result in results:
        buf = StringIO(result.stdout.decode('utf-8'))
        paf = parse_paf(buf)
        buf.close()
        print(paf)

if __name__ == "__main__":
    argv0 = argv[0]
    exit(main(argv[1:]))
