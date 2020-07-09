#!/usr/bin/env python3
"""
script to compare simulated breakpoints to those predicted by the algorithm
"""

import os
import sys
import json
import numpy as np
import subprocess

from sys  import argv, exit
from io import StringIO
from glob import glob
from itertools import chain
from collections import Counter

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

sys.path.insert(0, os.path.abspath('.')) # gross hack
from pangraph.utils import parse_paf, breakpoint

argv0 = None

def mode(items):
    return Counter(items).most_common(1)[0][0]

def rm_prefix(s, pfx):
    if s.startswith(pfx):
        return s[len(pfx):]
    return s

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
# [graph num]{block hash} -> [Match]
# key is the predicted blocks. value is the ancestral
class Matches():
    def __init__(self, pangraph, dir):
        with open(f"{dir}/ancestral.json") as fd:
            self.ancestral = json.load(fd)

        self.mu, self.beta = rm_prefix(pangraph, f"{dir}/").split(".")[:2]
        with open(pangraph) as fd:
            self.pangraph = json.load(fd)["tree"]["graph"]

        init = lambda g: {b["id"]:[] for b in g["blocks"]}
        self.graphs = [init(g) for g in self.pangraph]

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
        return [m[-1].pos[1] - m[0].pos[0] for g in self.graphs for m in g.values() if len(m) > 0]

    def coverage(self):
        found, hidden = [], []
        for n, g in enumerate(self.graphs):
            blocks = {b["id"]:b for b in self.pangraph[n]["blocks"]}
            for b, ms in g.items():
                if len(ms) == 0:
                    # TODO: log this in a better way?
                    hidden.append(None)
                    continue

                isolates = set([int(rm_prefix(k.split('?')[0],"isolate_")) for k in blocks[b]["muts"].keys()])
                # "found" ancestral blocks
                for m in [ms[0], ms[-1]]:
                    anc, blk = m.id
                    anc  = self.ancestral[str(anc)]["geneology"][blk]["present"]
                    date = [a["date"] for iso in isolates for a in anc[str(iso)]]
                    found.extend(date)
                # "hidden" ancestral blocks
                if len(ms) > 1:
                    for m in ms[1:-1]:
                        anc, blk = m.id
                        anc  = self.ancestral[str(anc)]["geneology"][blk]["present"]
                        date = [ a["date"] for iso in isolates for a in anc[str(iso)] ]
                        hidden.extend(date)

        return found, hidden

    def accuracy(self):
        return list( chain.from_iterable((m[0].diff[0], m[-1].diff[1]) for g in self.graphs for m in g.values() if len(m) > 0) )

def to_fastas(path):
    root = '/'.join(path.split('/')[:-1])
    with open(path, 'r') as ifd:
        pangraph = json.load(ifd)["tree"]["graph"]
        for n, graph in enumerate(pangraph):
            s  = [SeqRecord(id=b["id"], seq=Seq(b["seq"])) for b in graph["blocks"]]
            fa = f"{root}/graph_{n:03d}"
            with open(fa, 'w') as ofd:
                SeqIO.write(s, ofd, "fasta")
            yield fa

def usage():
    print(f"usage: {argv0} [directory]", file=sys.stderr)
    return 1

def main(args):
    if len(args) != 1 or not os.path.isdir(args[0]):
        exit(usage())
    dir      = args[0].rstrip("/")
    graphs   = glob(f"{dir}/*.*.pangraph.json")
    anc_blks = f"{dir}/ancestral.fa"

    if len(graphs) < 1:
        exit(f"directory {dir}: missing aligned pangraphs")
    if not os.path.exists(anc_blks):
        exit(f"directory {dir}: missing ancestral block sequences")

    all_matches = {}
    for n, g in enumerate(graphs):
        matches = Matches(g, dir)
        for i, sg in enumerate(to_fastas(g)):
            cmd = subprocess.run(
                ["minimap2", "-x", "asm5", str(anc_blks), sg],
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
                matches.put(i, hit['qry']['name'],
                    Match(*anc, hit['qry']['start'], hit['qry']['end'], d_beg, d_end, score)
                )

        matches.sort()
        all_matches[(matches.mu,matches.beta)] = matches

    def stats(matches):
        found, hidden = matches.coverage()
        return {'found':list(found), 'hidden':list(hidden), 'accuracy':list(matches.accuracy()), 'length':list(matches.length())}

    proc_matches = {f"({k[0]}, {k[1]})":stats(match) for k, match in all_matches.items()}
    with open(f"{dir}/algo_stats.json", "w") as fd:
        json.dump(proc_matches, fd)

    return 0

if __name__ == "__main__":
    argv0 = argv[0]
    exit(main(argv[1:]))
