#!/usr/bin/env python3
"""
script to measure the alignability of sequence around minimap2 breakpoints
"""

import os
import sys
import json

from glob import glob
from sys  import exit, argv
from Bio  import SeqIO

import numpy as np
import matplotlib.pylab as plt

import pyopa
from pyopa import Sequence as DNA
from pyopa import align_double as align

sys.path.insert(0, os.path.abspath('.')) # gross hack
from pangraph.tree  import Tree
from pangraph.utils import parse_paf

# ------------------------------------------------------------------------
# globals

EXTEND = 5000

# ------------------------------------------------------------------------
# helpers

def log(*msg):
    print(*msg, file=sys.stderr)

def seq_dict(path):
    with open(path, 'r') as rdr:
        d = {s.name:s.seq for s in SeqIO.parse(rdr, 'fasta')}
    return d

default = pyopa.load_default_environments()
env     = pyopa.generate_env(default['log_pam1'], 250)
def align_score(s1, s2):
    dna = lambda s: DNA(str(s))
    r   = align(dna(s1), dna(s2), env, False, True, True)
    return r[0]/max(len(s1), len(s2))

def count_blocks_in(hits, hit, delta):
    lb    = max(0, hit['start']-delta)
    ub    = min(hit['len'], hit['end']+delta)
    name  = hit['name']
    start = hit['start']
    end   = hit['end']

    num = 0
    for h in hits:
        if h['qry']['name'] == name:
            if lb <= h['qry']['start'] and h['qry']['start'] < start:
                num += 1
            elif end < h['qry']['start'] and h['qry']['start'] <= ub:
                num += 1
    return num

def draw_figure(ticks, left, right):
    Dl = np.diag(1/left[:,0])@left
    Dr = np.diag(1/right[:,0])@right
    fig, (ax1, ax2) = plt.subplots(1, 2)

    ax1.plot(ticks, Dl.T, color='r', alpha=.01)
    ax1.set_xlabel("extension (bp)")
    ax1.set_ylabel("alignment score (normalized)")
    ax1.set_ylim([0, 2])
    ax1.set_title("left junction")

    ax2.plot(ticks, Dr.T, color='b', alpha=.01)
    ax2.set_xlabel("extension (bp)")
    ax2.set_ylabel("alignment score (normalized)")
    ax2.set_ylim([0, 2])
    ax2.set_title("right junction")

    fig.savefig("figs/alignment_score_extend_past_junction.png", bbox_inches='tight')

# ------------------------------------------------------------------------
# main point of entry

def main(args):
    ticks     = np.arange(0, EXTEND, 250)
    ljunctions, rjunctions = [], []
    nblks_good = []
    nblks_all  = []
    for d in args:
        print(f"analyzing {d}")
        with open(f"{d}/guide.json") as fd:
            T = Tree.from_json(fd)

        path = lambda name, ext: f"{d}/tmp/{name}.{ext}"
        l_iv = lambda h: (h['start'], min(h['start']+1000, h['len']))
        r_iv = lambda h: (max(h['end']-1000, 0), h['end'])

        for n in T.preterminals():
            seq = {}
            for c in n.child:
                seq.update(seq_dict(path(c.name, "fa")))
            with open(path(n.name, "paf")) as fh:
                hits = parse_paf(fh)

            # TODO: use multiprocessing pool here...
            for hit in hits:
                qs, rs = seq[hit['qry']['name']], seq[hit['ref']['name']]
                ql, rl = l_iv(hit['qry']), l_iv(hit['ref'])
                qr, rr = r_iv(hit['qry']), r_iv(hit['ref'])
                lscore = lambda qb, rb: align_score(qs[max(qb,0):ql[1]], rs[max(rb,0):rl[1]])
                rscore = lambda qe, re: align_score(qs[qr[0]:min(qe, len(qs))], rs[rr[0]:min(re,len(rs))])

                ljunctions.append([lscore(ql[0]-tick, rl[0]-tick) for tick in ticks])
                rjunctions.append([rscore(qr[1]+tick, rr[1]+tick) for tick in ticks])

                if ljunctions[-1][-1] >= .9*ljunctions[-1][0] or rjunctions[-1][-1] >= .9*rjunctions[-1][0]:
                    # count number of blocks in extended zone
                    nblks_good.append(count_blocks_in(hits, hit['qry'], EXTEND))
                    nblks_good.append(count_blocks_in(hits, hit['ref'], EXTEND))
                    # for manual inspection
                    log(d, n.name)
                    log(json.dumps(hit))

                nblks_all.append(count_blocks_in(hits, hit['qry'], EXTEND))
                nblks_all.append(count_blocks_in(hits, hit['ref'], EXTEND))

    ljunctions, rjunctions = np.array(ljunctions), np.array(rjunctions)
    nblks_good, nblks_all  = np.array(nblks_good), np.array(nblks_all)
    np.savez("data/junctions.npz",
            left=ljunctions, right=rjunctions, ticks=ticks,
            blks_good=nblks_good, blks_all=nblks_all)
    draw_figure(ticks, ljunctions, rjunctions)

    return 0

if __name__ == "__main__":
    exit(main(argv[1:]))
