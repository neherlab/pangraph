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
from itertools import chain

from Bio  import SeqIO

sys.path.insert(0, os.path.abspath('.')) # gross hack
from pangraph.utils import parse_paf, breakpoint

argv0 = None

# simple struct
class Match():
    def __init__(self, ancestor, block, p_beg, p_end, d_beg, d_end, score):
        self.id    = (ancestor, block)
        self.pos   = (p_beg, p_end)
        self.diff  = (d_beg, d_end)
        self.score = score

    def __str__(self):
        return f"{{id: {self.id}; pos: {self.pos}; diff: {self.diff}; score: {self.score:.4f}}}"
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

    def put(self, n, blk, match):
        self.graphs[n][blk].append(match)

    def sort(self):
        for g in self.graphs:
            for key, val in g.items():
                g[key] = sorted(val, key=lambda m: m.pos[0])

    def length(self):
        return np.array([m[-1].pos[1] - m[0].pos[0] for g in self.graphs for m in g.values() if len(m) > 0])

    def coverage(self):
        return np.array([sum(len(ms)==1 for ms in g.values())/len(g) for g in self.graphs])

    def accuracy(self):
        return np.array(list(chain.from_iterable((m[0].diff[0], m[-1].diff[1]) for g in self.graphs for m in g.values() if len(m) > 0)))

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
    for n, g in enumerate(graphs):
        cmd = subprocess.run(
            ["minimap2", "-x", "asm5", str(anc_blks), str(g)],
            capture_output=True
        )
        buf  = StringIO(cmd.stdout.decode('utf-8'))
        hits = parse_paf(buf)
        buf.close()

        for hit in hits:
            anc   = [int(n) for n in hit['ref']['name'].split('_')[1:]]
            d_beg = hit['qry']['start']-hit['ref']['start']
            d_end = (hit['qry']['len']-hit['qry']['end']) - (hit['ref']['len']-hit['ref']['end'])
            score = hit['aligned_bases'] / hit['aligned_length']
            matches.put(n, hit['qry']['name'],
                Match(*anc, hit['qry']['start'], hit['qry']['end'], d_beg, d_end, score)
            )

    matches.sort()

    with open(f"{dir}/algo_stats.npz", "wb") as fd:
        np.savez(fd, coverage=matches.coverage(), accuracy=matches.accuracy(), length=matches.length())

    return 0

if __name__ == "__main__":
    argv0 = argv[0]
    exit(main(argv[1:]))
