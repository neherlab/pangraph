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

from pathos.multiprocessing import ProcessingPool as Pool

import numpy as np
import matplotlib.pylab as plt

from seqanpy import align_global as align

sys.path.insert(0, os.path.abspath('.')) # gross hack
from pangraph.tree  import Tree
from pangraph.utils import parse_paf, Strand, rev_cmpl

# ------------------------------------------------------------------------
# globals

SEQLEN = 1000
EXTEND = 5000

# ------------------------------------------------------------------------
# helpers

def log(*msg):
    print(*msg, file=sys.stderr)

def seq_dict(path):
    with open(path, 'r') as rdr:
        d = {s.name:s.seq for s in SeqIO.parse(rdr, 'fasta')}
    return d

default_opts = {
    'band'           : +1000,
    'score_match'    : +3,
    'score_mismatch' : -3,
    'score_gapopen'  : -5,
    'score_gapext'   : -1,
}
def do_align(s1, s2, aln_opts=default_opts):
    r   = align(s1, s2, **aln_opts)
    return {'score':max(r[0], 0)/max(len(s1), len(s2)), 'align' : r[1:]}

def align_score(s1, s2, aln_opts=default_opts):
    return do_align(s1, s2, aln_opts)['score']

def align_seqs(seqs, *ticks):
    class Aln:
        def __init__(self, qry, ref):
            aln = do_align(qry, ref)
            self.score = aln['score']
            self.seqs  = aln['align']
    class Hit:
        def __init__(self, seqs, delta):
            self.left  = Aln(seqs.left.qry(delta), seqs.left.ref(delta))
            self.right = Aln(seqs.right.qry(delta), seqs.right.ref(delta))

    def do(tick):
        return Hit(seqs, tick)

    if len(ticks) >= 5:
        with Pool(10) as jobs:
            return jobs.map(do, ticks)
    else:
        return [Hit(seqs, tick) for tick in ticks]

def write(wtr, *seqs):
    if len(seqs) == 0:
        return False

    for i, seq in enumerate(seqs):
        wtr.write(f">isolate_{i:02d}\n")
        wtr.write(f"{seq}\n")

    return True

def slice_seqs(qry, ref, hit):
    class Aln:
        def __init__(self):
            self.qry = None
            self.ref = None
    class Seqs:
        def __init__(self):
            self.left  = Aln()
            self.right = Aln()

    seqs = Seqs()

    qbeg, rbeg = hit['qry']['start'], hit['ref']['start']
    qend, rend = hit['qry']['end'],   hit['ref']['end']
    qlen, rlen = hit['qry']['len'],   hit['ref']['len']
    if hit['orientation'] == Strand.Plus:
        seqs.left.ref = lambda d: ref[max(rbeg-d, 0):min(rbeg+SEQLEN, qlen)]
        seqs.left.qry = lambda d: qry[max(qbeg-d, 0):min(qbeg+SEQLEN, qlen)]

        seqs.right.ref = lambda d: ref[max(rend-SEQLEN, 0):min(rend+d, rlen)]
        seqs.right.qry = lambda d: qry[max(qend-SEQLEN, 0):min(qend+d, qlen)]
    else:
        seqs.left.ref = lambda d: ref[max(rbeg-d, 0):min(rbeg+SEQLEN, qlen)]
        seqs.left.qry = lambda d: rev_cmpl(qry[max(qend-SEQLEN, 0):min(qend+d, qlen)])

        seqs.right.ref = lambda d: ref[max(rend-SEQLEN, 0):min(rend+d, rlen)]
        seqs.right.qry = lambda d: rev_cmpl(qry[max(qbeg-d, 0):min(qbeg+SEQLEN, qlen)])

    return seqs

def count_blocks_in(hits, hit, delta):
    lb    = max(0, hit['start']-delta)
    ub    = min(hit['len'], hit['end']+delta)
    name  = hit['name']
    start, end = hit['start'], hit['end']

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
    ax1.set_title("left junction")
    ax1.set_ylim([0, 2])

    ax2.plot(ticks, Dr.T, color='b', alpha=.01)
    ax2.set_xlabel("extension (bp)")
    ax2.set_ylabel("alignment score (normalized)")
    ax2.set_title("right junction")
    ax2.set_ylim([0, 2])

    fig.savefig("figs/alignment_score_extend_past_junction.png")

# ------------------------------------------------------------------------
# main point of entry

def main_align(args):
    ticks = np.arange(0, EXTEND, 250)
    ljunctions, rjunctions = [], []
    nblks_good, nblks_all  = [], []
    for d in args:
        print(f"analyzing {d}")
        with open(f"{d}/guide.json") as fd:
            T = Tree.from_json(fd)

        path = lambda name, ext: f"{d}/tmp/{name}.{ext}"
        for n in T.preterminals():
            seq = {}
            for c in n.child:
                seq.update(seq_dict(path(c.name, "fa")))

            with open(path(n.name, "paf")) as fh:
                hits = parse_paf(fh)

            for hit in hits:
                s   = slice_seqs(seq[hit['qry']['name']], seq[hit['ref']['name']], hit)
                aln = align_seqs(s, *ticks)

                ljunctions.append([a.left.score for a in aln])
                rjunctions.append([a.right.score for a in aln])

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

def main_scan(args):
    bands = list(range(0, 1000, 10))
    ngood_bps, nbad_bps = np.zeros(len(bands)), np.zeros(len(bands))
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

            for hit in hits:
                qs, rs = seq[hit['qry']['name']], seq[hit['ref']['name']]
                ql, rl = l_iv(hit['qry']), l_iv(hit['ref'])
                qr, rr = r_iv(hit['qry']), r_iv(hit['ref'])
                lscore = lambda qb, rb, opts: align_score(qs[max(qb,0):ql[1]], rs[max(rb,0):rl[1]], opts)
                rscore = lambda qe, re, opts: align_score(qs[qr[0]:min(qe, len(qs))], rs[rr[0]:min(re,len(rs))], opts)

                def do(b):
                    good, bad = 0, 0
                    opt = {
                        'band'           :  int(b),
                        'score_match'    : +3,
                        'score_mismatch' : -3,
                        'score_gapopen'  : -5,
                        'score_gapext'   : -1,
                    }
                    if lscore(ql[0]-EXTEND, rl[0]-EXTEND, opt) >= .75*lscore(ql[0], rl[0], opt):
                        bad += 1
                    else:
                        good += 1

                    if rscore(qr[1]+EXTEND, rr[1]+EXTEND, opt) >= .75*rscore(qr[1], rr[1], opt):
                        bad += 1
                    else:
                        good += 1
                    return (bad, good)

                with Pool(10) as jobs:
                    bad, good = zip(*jobs.map(do, bands))
                    nbad_bps  += np.array(bad)
                    ngood_bps += np.array(good)

    nbps = nbad_bps + ngood_bps
    fig  = plt.figure()
    ax   = plt.subplot(111)
    ax.plot(bands, nbad_bps/nbps)
    ax.set_xlabel("bandwidth")
    ax.set_ylabel("fraction of 'uniform' extension junctions")
    fig.savefig('figs/scan_breakpoints_bandwidth.png')

def main_gather(args):
    ROOT = "align"
    HI   = f"{ROOT}/1.5x"
    EQ   = f"{ROOT}/1.0x"
    LO   = f"{ROOT}/0.5x"

    NL = {'lo':0, 'hi':0, 'eq':0}
    NR = {'lo':0, 'hi':0, 'eq':0}
    if not os.path.exists(ROOT):
        os.mkdir(ROOT)
        os.mkdir(HI)
        os.mkdir(EQ)
        os.mkdir(LO)

    def proc(normal, extend, hit, prefix, N):
        if len(extend.seqs[0]) < EXTEND:
            return N

        if extend.score >= 1.5*normal.score:
            base = f"{HI}/{prefix}_{N['hi']}"
            indx = 'hi'
        elif extend.score >= 0.9*normal.score and extend.score <= 1.1*normal.score:
            base = f"{EQ}/{prefix}_{N['eq']}"
            indx = 'eq'
        elif extend.score <= 0.5*normal.score:
            base = f"{LO}/{prefix}_{N['lo']}"
            indx = 'lo'
        else:
            return N

        print(f"--> outputting {base}")
        with open(f"{base}.fa", 'w') as wtr:
            if write(wtr, *extend.seqs):
                N[indx] += 1
                with open(f"{base}.json", 'w') as wtr:
                    json.dump(hit, wtr, indent=4)
        return N

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

            for hit in hits:
                if hit['aligned_length'] < 1000:
                    continue

                s   = slice_seqs(seq[hit['qry']['name']], seq[hit['ref']['name']], hit)
                aln = align_seqs(s, 0, EXTEND)
                NL  = proc(aln[0].left,  aln[1].left,  hit, 'left',  NL)
                NR  = proc(aln[0].right, aln[1].right, hit, 'right', NR)

def main(args):
    if args[0] == "align":
        main_align(args[1:])
    elif args[0] == "scan":
        main_scan(args[1:])
    elif args[0] == "gather":
        main_gather(args[1:])

    return 0

if __name__ == "__main__":
    exit(main(argv[1:]))
