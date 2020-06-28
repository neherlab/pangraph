#!/usr/bin/env python3
"""
script to compare simulated breakpoints to those predicted by the algorithm
"""

import os
import sys
import numpy as np
import subprocess

from sys  import argv, exit
from io import StringIO
from glob import glob

from Bio  import SeqIO

sys.path.insert(0, os.path.abspath('.')) # gross hack
from pangraph.utils import parse_paf

argv0 = None

# simple struct
class Match():
    def __init__(self, ancestor, block, d_beg, d_end, score):
        self.id    = (ancestor, block)
        self.diff  = (d_beg, d_end)
        self.score = score

    def __str__(self):
        return f"{{id: {self.id};\tdiff: {self.diff};\tscore: {self.score:.4f}}}"
    def __repr__(self):
        return self.__str__()

# associative data structure
# [graph num]{block hash} -> Match
# key is the predicted blocks. value is the ancestral
class Matches():
    def __init__(self, graphs):
        self.graphs = [self.init(g) for g in graphs]

    def init(self, graph):
        return {s.id:[] for s in SeqIO.parse(graph, 'fasta')}

    def __str__(self):
        return "\n".join(str(g) for g in self.graphs)

    def __repr__(self):
        return self.__str__()

    def put(self, id, blk, match):
        self.graphs[id][blk].append(match)

    def coverage(self):
        return np.array([sum(len(ms)==1 for ms in g.values())/len(g) for g in self.graphs])

    def bp_accuracy(self):
        return np.array([m.diff[0] for g in self.graphs for ms in g.values() for m in ms])

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

    matches = Matches(graphs)
    for i, g in enumerate(graphs):
        cmd = subprocess.run(
            ["minimap2", "-x", "asm5", str(anc_blks), str(g)],
            capture_output=True
        )
        buf  = StringIO(cmd.stdout.decode('utf-8'))
        hits = parse_paf(buf)
        buf.close()

        for hit in hits:
            anc   = [int(n) for n in hit['ref']['name'].split('_')[1:]]
            d_beg = hit['qry']['start'] - hit['ref']['start']
            d_end = hit['qry']['end'] - hit['ref']['end']
            score = hit['aligned_bases'] / hit['aligned_length']
            matches.put(i, hit['qry']['name'],
                Match(*anc, d_beg, d_end, score)
            )

        print(matches.coverage())
        print(matches.bp_accuracy())

if __name__ == "__main__":
    argv0 = argv[0]
    exit(main(argv[1:]))
