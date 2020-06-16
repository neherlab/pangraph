import os

import numpy as np
import numpy.random as rng

from Bio import SeqIO
from portion import Interval, closedopen, to_data

from .graph import Graph
from .utils import asrecord, cat

# ------------------------------------------------------------------------
# rng sampling

def random_seq(L, alphabet=None, barcode=None):
    if alphabet is None:
        alphabet = np.frombuffer(b"ACGT", dtype=np.int8)

    idx = rng.randint(len(alphabet), size=L)
    seq = alphabet[idx]
    if barcode is None:
        return seq

    return seq, barcode*np.ones(L), np.arange(L)

def indel(seq, other=None, std=0, L0=0, N=0):
    if len(seq[0]) != len(seq[1]):
        print("unequal sequence lengths")
        import ipdb; ipdb.set_trace()

    L = len(seq[0])
    i = rng.randint(L)
    if other is not None:
        k = rng.randint(len(other[0]))

    dl  = int(rng.randn(1)*std + (L0-L))

    def interval(seq, i, j):
        if i < j:
            return seq[i:j]
        else:
            return cat(seq[i:], seq[:j])

    def do(s, o, isseq=True):
        if -L < dl < 0:
            j = (i - dl) % L
            if i < j:
                s = cat(s[:i], s[j:])
            else:
                s = s[j:i]
        elif 0 < dl < L:
            if o is None:
                ins = random_seq(dl) if isseq else N*np.ones(dl, dtype=np.int8)
            elif dl < len(o):
                j = (k + dl) % len(o)
                ins = interval(o, k, j)
            s = cat(s[:i], ins, s[i:])

        return s

    if other is None:
        seq = tuple(do(s, None, i==0) for i, s in enumerate(seq))
    else:
        seq = tuple(do(s, o, i==0) for i, (s, o) in enumerate(zip(seq, other)))

    return seq

def transpose(seq, avg, std):
    L  = len(seq)
    dl = int(rng.randn(1)*std + avg)

    if 0 < dl < L:
        i = rng.randint(L)
        k = rng.randint(len(seq[0]))

        def do(s):
            j = (i+dl) % L

            xfer = []
            if i < j:
                xfer = s[i:j]
                s    = cat(s[:i], s[j:])
            else:
                xfer = cat(s[i:], s[:j])
                s    = s[j:i]

            s = cat(s[:k], xfer, s[k:])
            return s

        seq = tuple(do(s) for s in seq)

    return seq

def invert(seq, avg, std):
    L  = len(seq)
    dl = int(rng.randn(1)*std + avg)
    if 0 < dl < L:
        i = rng.randint(L)
        j = (i+dl) % L

        def do(s):
            xfer = ""
            if i < j:
                s = cat(s[:i], s[i:j][::-1], s[j:])
            else:
                xfer = cat(s[i:], s[:j])
                xfer = xfer[::-1]
                s    = cat(xfer[:j], s[j:i], xfer[j:])

            return s

        seq = tuple(do(s) for s in seq)

    return seq

def mutate(seq, mu, alphabet=None):
    if alphabet is None:
        alphabet = np.array([ord(c) for c in ['A', 'C', 'G', 'T']])

    # Unpack
    seq, bc, anc = seq

    nm  = rng.poisson(lam=mu*len(seq), size=1)
    idx = rng.choice(len(seq), size=nm, replace=False)
    mut = rng.choice(len(alphabet), size=nm, replace=True)

    for i, I in enumerate(idx):
        while alphabet[mut[i]] == seq[I]:
            mut[i] = rng.choice(len(alphabet), 1)

        seq[I] = alphabet[mut[i]]

    # Repack
    return tuple((seq, bc, anc))

# ------------------------------------------------------------------------
# an associative data structure whose keys are 1d intervals

class IntervalMap(object):
    def __init__(self, low, hi):
        self.ival = []
        self.data = []
        self.boundary = [low, hi]

    def put(self, key, val, top_level=True):
        K, i = key, 0
        if not isinstance(val, set):
            val = set([val])

        while i < len(self.ival) and not K.empty:
            ival = self.ival[i]
            if not ival.overlaps(K):
                i += 1
                continue

            # terminate
            if ival == K:
                K = K - ival
                self.data[i].update(val)
                break

            # terminate
            elif ival.contains(K):
                Is = [Interval(iv) for iv in ival - K]
                Is.append(K)

                Ds = [self.data[i].copy() for _ in range(len(Is))]
                Ds[-1].update(val)

                idx = list(range(len(Is)))
                idx.sort(key=lambda x: Is[x].lower)
                Is, Ds = [Is[i] for i in idx], [Ds[i] for i in idx]

                self.ival = self.ival[:i] + Is + self.ival[i+1:]
                self.data = self.data[:i] + Ds + self.data[i+1:]

                K = K - ival
                break

            # add value to ival entry, remove from K, continue
            elif K.contains(ival):
                K = K - ival
                self.data[i].update(val)
                i += 1

            # partial hit on both
            else:
                tmp, K = K, K - ival
                ov = tmp.intersection(ival)
                Is = [ival-tmp]
                Ds = [self.data[i]]
                Ds[0].update(val)

                if not ov.empty:
                    Is.append(Interval(ov))
                    Ds.append(set(Ds[0]))

                    idx = list(range(len(Is)))
                    idx.sort(key=lambda x: Is[x].lower)
                    Is = [Is[x] for x in idx]
                    Ds = [Ds[x] for x in idx]

                self.ival = self.ival[:i] + Is + self.ival[i+1:]
                self.data = self.data[:i] + Ds + self.data[i+1:]

                i += len(Is)

        if not K.empty:
            for k in K:
                I = 0
                while I < len(self.ival) and self.ival[I].lower < k.lower:
                    I += 1
                self.ival.insert(I, Interval(k))
                self.data.insert(I, val)

        # find intervals that are fully translated by L
        # rotate them back down into the primary interval
        ok = True
        for i, I in enumerate(self.ival):
            if I.lower >= self.boundary[1]:
                ok = False
                break

        # NOTE: Assumes ivals are sorted by lower bounds!!
        if not ok:
            self.ival, ival_add = self.ival[:i], self.ival[i:]
            self.data, data_add = self.data[:i], self.data[i:]

            ival_add = [closedopen(iv.lower % self.boundary[1], iv.upper % self.boundary[1]) for iv in ival_add]
            for iv, d in zip(ival_add, data_add):
                self.put(iv, d, False)

        # find intervals that span the gap between primary & secondary intervals
        # XXX: not sure why this code is required...
        ok = True
        for i, I in enumerate(self.ival):
            if I.upper > self.boundary[1]:
                ok = False
                break

        if not ok and top_level:
            for iv, d in zip(self.ival[:i], self.data[:i]):
                I = closedopen(iv.lower+self.boundary[1], iv.upper+self.boundary[1])
                self.put(I, d, False)

        return

    def get(self, key):
        data = []
        for i, ival in enumerate(self.ival):
            if ival.overlaps(key):
                data.append((ival, self.data[i]))

        return data

    def serialize(self):
        def unpack(datum):
            keys, vals = [], []
            for d in datum:
                keys.append(int(d[0]))
                vals.append(tuple(int(e) for e in d[1]))
            return dict(zip(keys,vals))

        def pack(intervals):
            return [tuple([int(e) for e in d[1:3]] for d in to_data(iv)) for iv in intervals]

        if len(self.ival) == 0:
            return None

        return {
            "interval" : pack(self.ival),
            "children" : [unpack(d) for d in self.data],
        }

# ------------------------------------------------------------------------
# main class used to generate data

class Population(object):
    def __init__(self, size, len, mu, rate_hgt, rate_indel, rate_transpose):
        # Population parameters
        self.N    = size
        self.L    = len
        self.mu   = mu

        # reassortment parameters
        # TODO: make adjustable
        self.seqmv_avg = self.L // 5
        self.seqmv_std = self.seqmv_avg // 2
        self.indel_std = self.L // 20

        # reassortment rates (per generation)
        # perform cumulative sum here
        self.rhgt       = rate_hgt
        self.rindel     = self.rhgt + rate_indel
        self.rtranspose = self.rindel + rate_transpose

        # store thunks for nicer function calls
        self.hgt       = lambda seq, donor: indel(seq, donor, self.indel_std, self.L, self.N)
        self.indel     = lambda seq: indel(seq, None, self.indel_std, self.L, self.N)
        self.transpose = lambda seq: transpose(seq, self.seqmv_avg, self.seqmv_std)
        self.mutate    = lambda seq: mutate(seq, self.mu)

        # Population data
        self.seq = [random_seq(self.L, barcode=n) for n in range(self.N)]
        self.anc = self.seq.copy()

    def evolve(self, T):
        for _ in range(T):
            rand   = rng.rand(self.N)
            parent = rng.choice(self.N, size=self.N)
            seq    = []
            for n in range(self.N):
                s = self.seq[parent[n]]
                # NOTE: this must be in the same order as the rates above
                if rand[n] < self.rhgt:
                    donor = self.seq[int(rng.choice(self.N, size=1))]
                    s = self.hgt(s, donor)
                elif rand[n] < self.rindel:
                    s = self.indel(s)
                elif rand[n] < self.rtranspose:
                    s = self.transpose(s)

                s = self.mutate(s)
                seq.append(s)

            self.seq = seq

    def write_fasta(self, wtr, prefix="isolate"):
        seq   = lambda s: "".join(chr(n) for n in s[0])
        name  = lambda n: f"{prefix}_{n:05d}"
        entry = [asrecord(seq(s), name(n)) for n, s in enumerate(self.seq)]

        return SeqIO.write(entry, wtr, "fasta")

    def ancestral_weights(self):
        W = np.zeros((N, L))
        for _, anc, pos in self.seq:
            for a, x in zip(anc, pos):
                if a < N:
                    W[int(a), int(x)] += 1

        return W

    def ancestral_blocks(self):
        blks = [IntervalMap(0, self.L) for _ in range(len(self.seq))]

        # ---------------------------------
        # internal funcs

        def And(x, y):
            return np.bitwise_and(x, y)

        def Put(A, anciv, C, curiv, l):
            if curiv[0] >= curiv[1]:
                curiv = (curiv[0], curiv[1]+l)
            vals = [(C, closedopen(curiv[0], curiv[1]))]

            for val in vals:
                if anciv[0] >= anciv[1]:
                    anciv = (anciv[0], anciv[1]+self.L)
                blks[int(A)].put(closedopen(anciv[0], anciv[1]+1), val)

        # ---------------------------------
        # body

        for n, (_, anc, pos) in enumerate(self.seq):
            delta     = np.empty(pos.shape)
            delta[0]  = np.abs(pos[0]  - pos[-1])
            delta[1:] = np.abs(pos[1:] - pos[:-1])

            cond  = And(And(delta != 1, delta != self.L-1), anc < self.N)
            bkpnt = np.where(cond)[0]

            if len(bkpnt) == 0:
                # TODO: Put check for difference in ancestral values here??
                Put(anc[0], (pos[0], pos[-1]), n, (0, len(pos)), len(anc))
                continue

            for i, bp in enumerate(bkpnt[:-1]):
                nbp = bkpnt[i+1]
                assert bp < nbp, "non-ordered breakpoints"
                if not np.equal(anc[bp:nbp], anc[bp]).all():
                    panic("inconsistent cut/boundary")
                else:
                    Put(anc[bp], (pos[bp], pos[nbp-1]), n, (bp, nbp), len(anc))

            bp, nbp = bkpnt[-1], bkpnt[0]
            Put(anc[bp], (pos[bp], pos[nbp-1]),  n, (bp, nbp), len(anc))

        # Now that we've found the ancestral block 'units', we go back and update the intervals
        for nanc, blk in enumerate(blks):
            for i in range(len(blk.ival)):
                anc_iv, data = blk.ival[i], blk.data[i]
                updated = set([])
                for item in data:
                    ncur, _ = item
                    anc_id, anc_pos = self.seq[ncur][1:3]
                    if len(anc_iv) > 1:
                        import ipdb; ipdb.set_trace()

                    # TODO: Make faster!
                    lbs = np.where(And(anc_pos == anc_iv.lower, anc_id==nanc))[0]
                    ubs = np.where(And(anc_pos == (anc_iv.upper-1)%self.L, anc_id==nanc))[0]
                    for lb, ub in zip(sorted(lbs), sorted(ubs)):
                        updated.add((ncur, (lb, ub+1)))

                blk.data[i] = updated

        # Return the sparse array
        sparse = {}
        for i, b in enumerate(blks):
            data = b.serialize()
            if data is not None:
                sparse[i] = data

        return sparse

# ------------------------------------------------------------------------
# Main point of entry

# def test_interval_map():
#     print("ROUND 1")
#     M = IntervalMap()
#     M.put(IV.closedopen(0, 10), 1)
#     M.put(IV.closedopen(2, 5), 2)
#     M.put(IV.closedopen(12, 15), 3)
#     M.put(IV.closedopen(30, 45), 3)
#     M.put(IV.closedopen(0, 100), 4)
#     print(f"--> Keys: {M.ival}")
#     print(f"--> Vals: {M.data}")
#     print(f"--> Return: {M.get(IV.closedopen(1, 5))}")

#     print("ROUND 2")
#     M = IntervalMap()
#     M.put(IV.closedopen(30, 45), 1)

#     x = IV.closedopen(70, 100)
#     y = IV.closedopen(0, 20)
#     z = x.union(y)
#     M.put(z, 2)

#     print(f"--> Keys: {M.ival}")
#     print(f"--> Vals: {M.data}")

#     M.put(IV.closedopen(15, 35), 3)

#     print(f"--> Keys: {M.ival}")
#     print(f"--> Vals: {M.data}")

# NWK_DIR = "data/synth/nwk"
# MTX_DIR = "data/synth/mtx"
# PKL_DIR = "data/synth/pkl"
# SEQ_DIR = "data/synth/seq"
# VIS_DIR = "data/synth/vis"
# BLK_DIR = "data/synth/blk"

# NITS = 10
# MUS  = np.linspace(0, 1e-6, 20)
# TS   = np.array([ int(N / n) for n in [7, 5, 3, 1, .5, .25] ], dtype=np.int)

# def mk_grid():
#     for M, mu in enumerate(MUS):
#         for I, T in enumerate(TS):
#             for it in range(NITS):
#                 fa   = f"{SEQ_DIR}/mu{M:02d}_T{I:02d}_{it:02d}.fa"
#                 bfa  = f"{BLK_DIR}/mu{M:02d}_T{I:02d}_{it:02d}.fa"
#                 mtx  = f"{MTX_DIR}/mu{M:02d}_T{I:02d}_{it:02d}.mtx"
#                 tree = f"{NWK_DIR}/mu{M:02d}_T{I:02d}_{it:02d}.nwk"

#                 print(f"Building random seqs")
#                 seqs = [random_seq(L, barcode=n) for n in range(N)]
#                 ancs = seqs
#                 for _ in range(T):
#                     seqs = evolve(seqs, mu)

#                 # W    = ancestral_weights(seqs)
#                 print(f"Collecting blocks")
#                 blks = ancestral_blocks(seqs)

#                 # Get distribution of block lengths
#                 lens, nums, blkseqs = [], [], []
#                 for n, blk in enumerate(blks):
#                     for i, iv in enumerate(blk.ival):
#                         lens.append(iv.upper-iv.lower)
#                         nums.append(len(blk.data[i]))
#                         for atom in iv:
#                             if atom.upper <= L0:
#                                 blkseqs.append([ancs[n][0][atom.lower:atom.upper]])
#                             else:
#                                 blkseqs.append([np.array(ancs[n][0][atom.lower:L0].tolist() + \
#                                                          ancs[n][0][0:atom.upper].tolist() )])
#                 lens = np.array(lens)

#                 print(f"Number of blocks: {len(lens)}")
#                 # print(f"--> Output {bfa}")
#                 with open(bfa, "w+") as out:
#                     tofasta(blkseqs, out, prefix="block")

#                 # xratio = np.sum(lens/np.sum(lens) * nums)
#                 # print(f"Compression Ratio: {xratio}")
#                 # cdfplot(lens)

#                 # Output individual current sequences
#                 # print(f"--> Output {fa}")
#                 with open(fa, "w+") as out:
#                     tofasta(seqs, out)

#                 # Output kmer guide tree
#                 mash(fa, mtx)
#                 # print(f"--> Output {mtx}")
#                 with open(mtx) as out:
#                     D, names = getmtx(out)

#                 sf  = squareform(D)
#                 Z   = linkage(sf)
#                 t   = to_tree(Z)
#                 nwk = getnwk(t, "", t.dist, names)

#                 # print(f"--> Output {tree}")
#                 with open(tree, 'w+') as out:
#                     out.write(f"{nwk}\n")

#     with open("data/synth/params.npz", "w+") as out:
#         np.savez(out, Mus=MUS, Ts=TS)

#     # # Run graph creation code
#     # f = "test.nwk"
#     # Ss = fai.Fasta(fa)

#     # g, ok = Graph.fromnwk(f, Ss) #, mu=0, beta=0)
#     # g.name = "test"
#     # g.tofasta()
#     # pickle.dump(g.todict(), open(f"{g.name}.pkl", "wb"))

#     # os.system(f"minimap2 -t 2 -x asm5 -D -c blocks.fasta data/graph/aln/test.fasta 1>test.paf 2>log")
#     # hits = parsepaf("test.paf")

# if __name__ == "__main__":
#     # --------------------
#     # Possible entry points
#     # test_interval_map()
#     # mk_grid()
#     # --------------------

#     import re
#     import matplotlib.pylab as plt
#     from py.utils import cdfplot

#     T  = N // 7
#     mu = 0

#     STOP = False
#     for nit in range(20):
#         if STOP:
#             break
#         print(f"Running iteration {nit:02d}")
#         seqs = [random_seq(L, barcode=n) for n in range(N)]
#         ancs = seqs.copy()
#         for _ in range(T):
#             seqs = evolve(seqs, mu)

#         print(f"Collecting blocks")
#         blks = ancestral_blocks(seqs)

#         blkseqs = {}
#         for n, blk in enumerate(blks):
#             for i, iv in enumerate(blk.ival):
#                 for j, atom in enumerate(iv):
#                     assert atom.lower < L0, "bad interval"
#                     if atom.upper < L0:
#                         blkseqs[(n, i, j)] = "".join(chr(c) for c in ancs[n][0][atom.lower:atom.upper])
#                     else:
#                         blkseqs[(n, i, j)] = "".join(chr(c) for c in (ancs[n][0][atom.lower:L0].tolist() + \
#                                                                       ancs[n][0][0:atom.upper].tolist()))

#         # Convert sequences to strings
#         for i, seq in enumerate(seqs):
#             seqs[i] = ("".join(chr(c) for c in seq[0]), seq[1], seq[2])

#         for (b, i, j), blk in blkseqs.items():
#             if len(blk) < 10:
#                 continue

#             locs = []
#             for n, seq in enumerate(seqs):
#                 locs.extend([(n, m.start()) for m in re.finditer(blk, seq[0])])

#             data = set()
#             for elt in blks[b].data[i]:
#                 if elt[1][0] < elt[1][1]:
#                     if abs(elt[1][1]-elt[1][0] - (len(blk)-1)) < 10:
#                         data.add((elt[0], elt[1][0]))
#                 else:
#                     if abs((len(seqs[elt[0]][0]) - elt[1][0]) + elt[1][1] - (len(blk)-1)) < 10:
#                         data.add((elt[0], elt[1][1]))

#             locs = set(locs)
#             if locs != data and len(blk) > 100:
#                 print(locs-data)
#                 STOP = True
#                 break
