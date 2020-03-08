import os
import pickle
import warnings
warnings.simplefilter("ignore", DeprecationWarning)

import numpy as np
import numpy.random as rng

from scipy.cluster.hierarchy import dendrogram, linkage, to_tree
from scipy.spatial.distance import squareform

import intervals as IV
import pyfaidx   as fai

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from py.utils import asstring, parsepaf
from py.graph import Graph

# ------------------------------------------------------------------------
# Class

class Population(object):

    def __init__(self, N, L, mu, dir, name):
        # Population parameters
        self.N    = N
        self.L    = L
        self.mu   = mu
        self.dir  = dir
        self.name = name

        # Distribution parameters
        self.lx = L // 5
        self.ls = self.lx // 2
        self.Ls = L0 // 20

        # Population data
        self.seqs = [random_seq(L, barcode=n) for n in range(N)]

    def evolve(self, T):
        for _ in range(T):
            pass

# ------------------------------------------------------------------------
# Globals

N  = int(100)
L  = int(1e5)

L0 = L
Ls = L0 // 20

lx = L // 5
ls = lx // 2

# ------------------------------------------------------------------------
# Useful structs

def atomic(iv):
    if isinstance(iv, IV.AtomicInterval):
        return iv
    return iv.to_atomic()

def tostr(x):
    return "".join(chr(e) for e in x)

class IntervalMap(object):
    def __init__(self):
        self.ival = []
        self.data  = []


    def put(self, key, val, top_level=True):
        K, i = key, 0
        if not isinstance(val, set):
            val = set([val])

        while i < len(self.ival) and not K.is_empty():
            ival = self.ival[i]
            if not ival.overlaps(K):
                i += 1
                continue

            if ival == K:
                K = K - ival
                self.data[i].update(val)
                break

            elif ival.contains(K):
                Is = [IV.Interval(iv) for iv in ival - K]
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

            elif K.contains(ival):
                K = K - ival
                # if len(K) == 1:
                self.data[i].update(val)
                i += 1
                # else:
                #     Ds = [set([val]) for _ in range(len(K))]
                #     Ks = [IV.Interval(k) for k in K[:-1]] + [ival]

                #     self.ival = self.ival[:i] + Ks + self.ival[i+1:]
                #     self.data = self.data[:i] + Ds + self.data[i+1:]
                #     i += len(Ks)

                #     K = IV.Interval(K[-1])
            else:
                tmp, K = K, K - ival
                ov = tmp.intersection(ival)
                Is = [ival-tmp]
                Ds = [self.data[i]]
                Ds[0].update(val)

                if not ov.is_empty():
                    Is.append(IV.Interval(ov))
                    Ds.append(set(Ds[0]))

                    idx = list(range(len(Is)))
                    idx.sort(key=lambda x: Is[x].lower)
                    Is = [Is[x] for x in idx]
                    Ds = [Ds[x] for x in idx]

                self.ival = self.ival[:i] + Is + self.ival[i+1:]
                self.data = self.data[:i] + Ds + self.data[i+1:]

                i += len(Is)

        if not K.is_empty():
            for k in K:
                I = 0
                while I < len(self.ival) and self.ival[I].lower < k.lower:
                    I += 1
                self.ival.insert(I, IV.Interval(k))
                self.data.insert(I, val)

        # Find intervals that are fully translated by L0
        ok = True
        for i, I in enumerate(self.ival):
            if I.lower >= L0:
                ok = False
                break

        # NOTE: Assumes ivals are sorted by lower bounds!!
        if not ok:
            self.ival, ival_add = self.ival[:i], self.ival[i:]
            self.data, data_add = self.data[:i], self.data[i:]

            ival_add = [IV.closedopen(iv.lower%L0, iv.upper%L0) for iv in ival_add]
            for iv, d in zip(ival_add, data_add):
                self.put(iv, d, False)

        ok = True
        for i, I in enumerate(self.ival):
            if I.upper > L0:
                ok = False
                break

        if not ok and top_level:
            # import ipdb; ipdb.set_trace()
            for iv, d in zip(self.ival[:i], self.data[:i]):
                I = IV.openclosed(iv.lower+L0, iv.upper+L0)
                self.put(I, d, False)
            # import ipdb; ipdb.set_trace()


        return None

    def get(self, key):
        data = []
        for i, ival in enumerate(self.ival):
            if ival.overlaps(key):
                data.append((ival, self.data[i]))

        return data

# ------------------------------------------------------------------------
# Functions

def mash(fasta, out):
    cmd = f"~/opt/mash-Linux64-v2.2/mash triangle {fasta} > {out}"
    os.system(cmd)

def getmtx(file):
    N = int(file.readline().strip())
    D, names = np.zeros((N, N)), []
    for i, line in enumerate(file):
        column = line.split()
        names.append(column[0])
        for j, data in enumerate(column[1:]):
            D[i, j] = float(data)
            D[j, i] = D[i, j]

    return D, names

def cat(*args):
    return np.concatenate(tuple(arg for arg in args))

def grab_interval(seq, i, j):
    if i < j:
        return seq[i:j]
    else:
        return cat(seq[i:], seq[:j])

def getnwk(node, newick, parentdist, leaf_names):
    if node.is_leaf():
        return "%s:%.8f%s" % (leaf_names[node.id], parentdist - node.dist, newick)
    else:
        if len(newick) > 0:
            newick = "):%.8f%s" % (parentdist - node.dist, newick)
        else:
            newick = ");"
        newick = getnwk(node.get_left(), newick, node.dist, leaf_names)
        newick = getnwk(node.get_right(), ",%s" % (newick), node.dist, leaf_names)
        newick = "(%s" % (newick)
        return newick

def random_seq(L, alphabet=None, barcode=None):
    if alphabet is None:
        alphabet = np.frombuffer(b"ACGT", dtype=np.int8)

    idx = rng.randint(len(alphabet), size=L)
    seq = alphabet[idx]
    if barcode is None:
        return seq

    return seq, barcode*np.ones(L), np.arange(L)

def random_indel(seq, other=None):
    if len(seq[0]) != len(seq[1]):
        print("unequal sequence lengths")
        import ipdb; ipdb.set_trace()

    L = len(seq[0])
    i = rng.randint(L)
    k = rng.randint(len(other[0]))

    mu = L0 - L
    dl = int(rng.randn(1)*Ls + mu)

    def do(s, o, isseq=True):
        if -L < dl < 0:
            j = (i - dl) % L
            if i < j:
                s = cat(s[:i], s[j:])
            else:
                s = s[j:i]
        elif 0 < dl < L:
            if o is None:
                ins = random_seq(dl) if isseq else N*np.ones(dl)
            elif dl < len(o):
                j = (k + dl) % len(o)
                ins = grab_interval(o, k, j)
            s = cat(s[:i], ins, s[i:])

        return s

    if other is None:
        seq = tuple(do(s, None, i==0) for i, s in enumerate(seq))
    else:
        seq = tuple(do(s, o, i==0) for i, (s, o) in enumerate(zip(seq, other)))

    return seq

def random_transposition(seq):
    L  = len(seq)
    dl = int(rng.randn(1)*ls + lx)

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

def random_invert(seq):
    L  = len(seq)
    dl = int(rng.randn(1)*ls + lx)
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
        alphabet = "ACGT"

    # Unpack
    seq, bc, anc = seq

    nm  = rng.poisson(lam=int(mu*len(seq)), size=1)
    idx = rng.choice(len(seq), size=nm, replace=False)
    mut = rng.choice(len(alphabet), size=nm, replace=False)

    for i, I in enumerate(idx):
        while alphabet[mut[i]] == seq[I]:
            mut[i] = rng.choice(len(alphabet), 1)

        seq[i] = alphabet[mut[i]]

    # Repack
    return tuple((seq, bc, anc))

r_transposition  = .00
r_indel = .00
r_hgt  = .10
def evolve(seqs, mu=0):
    nseqs   = []
    parent  = rng.choice(len(seqs), size=len(seqs))
    for i in range(len(seqs)):
        nseq = seqs[parent[i]]
        rand = rng.rand(1)
        if rand < r_transposition:
            nseq = random_transposition(nseq)
        elif rand < r_indel + r_transposition:
            nseq = random_indel(nseq)
        elif rand < r_indel + r_transposition + r_hgt:
            oseq = seqs[int(rng.choice(len(seqs), size=1))]
            nseq = random_indel(nseq, oseq)

        if mu > 0:
            nseq = mutate(nseq, mu)
        nseqs.append(nseq)

    return nseqs

def tofasta(seqs, wtr, prefix="isolate"):
    recs = []
    for n, seq in enumerate(seqs):
        seq  = "".join(chr(n) for n in seq[0])
        name = f"{prefix}_{n:04d}"
        recs.append(SeqRecord(Seq(seq), id=name, name=name))

    return SeqIO.write(recs, wtr, "fasta")

def ancestral_weights(seqs):
    W = np.zeros((N, L))
    for _, anc, pos in seqs:
        for a, x in zip(anc, pos):
            if a < N:
                W[int(a), int(x)] += 1

    return W

def ancestral_blocks(seqs):
    # NOTE: Important -- We assume the ancestral plasmid length is L0. If this changes,
    #                    this code is wrong as the modulo arithmetic needs to change.
    blks = [IntervalMap() for _ in range(len(seqs))]

    # ---------------------------------
    # Internal Functions

    def And(x, y):
        return np.bitwise_and(x, y)

    def Put(A, anciv, C, curiv, l):
        if curiv[0] >= curiv[1]:
            curiv = (curiv[0], curiv[1]+l)
        vals = [(C, IV.closedopen(curiv[0], curiv[1]))]

        # if int(A) == 44:
        #     print(f"{anciv}->{blks[int(A)].ival}")
        #     print(f"--> iso: {C}, val: {vals}")
        #     import ipdb; ipdb.set_trace()

        for val in vals:
            if anciv[0] >= anciv[1]:
                anciv = (anciv[0], anciv[1]+L0)
            blks[int(A)].put(IV.closedopen(anciv[0], anciv[1]+1), val)

        # if int(A) == 44:
        #     print(f"RESULT: {blks[int(A)].ival}")
        #     import ipdb; ipdb.set_trace()
    # ---------------------------------

    for n, (_, anc, pos) in enumerate(seqs):
        delta     = np.empty(pos.shape)
        delta[0]  = np.abs(pos[0]  - pos[-1])
        delta[1:] = np.abs(pos[1:] - pos[:-1])

        cond  = And(And(delta != 1, delta != L0-1), anc < N)
        bkpnt = np.where(cond)[0]

        if len(bkpnt) == 0:
            # TODO: Put check for difference in ancestral values here??
            Put(anc[0], (pos[0], pos[-1]), n, (0, len(pos)), len(anc))
            continue

        for i, bp in enumerate(bkpnt[:-1]):
            nbp = bkpnt[i+1]
            assert bp < nbp, "non-ordered breakpoints"
            if not np.equal(anc[bp:nbp], anc[bp]).all():
                print("PANIC: Inconsistent cut/boundary")
                import ipdb; ipdb.set_trace()
            else:
                Put(anc[bp], (pos[bp], pos[nbp-1]), n, (bp, nbp), len(anc))

        bp, nbp = bkpnt[-1], bkpnt[0]
        Put(anc[bp], (pos[bp], pos[nbp-1]),  n, (bp, nbp), len(anc))

    # Now that we've found the ancestral block 'units', we go back and update the intervals
    for nanc, blk in enumerate(blks):
        for i in range(len(blk.ival)):
            # if nanc == 11:
            #     import ipdb; ipdb.set_trace()

            anc_iv, data = blk.ival[i], blk.data[i]
            updated = set([])
            for item in data:
                ncur, _ = item
                anc_id, anc_pos = seqs[ncur][1:3]
                if len(anc_iv) > 1:
                    import ipdb; ipdb.set_trace()

                # TODO: Make faster!
                lbs = np.where(And(anc_pos == anc_iv.lower, anc_id==nanc))[0]
                ubs = np.where(And(anc_pos == (anc_iv.upper-1)%L0, anc_id==nanc))[0]
                for lb, ub in zip(sorted(lbs), sorted(ubs)):
                    updated.add((ncur, (lb, ub)))

            blk.data[i] = updated

    return blks

# ------------------------------------------------------------------------
# Main point of entry

def test_interval_map():
    print("ROUND 1")
    M = IntervalMap()
    M.put(IV.closedopen(0, 10), 1)
    M.put(IV.closedopen(2, 5), 2)
    M.put(IV.closedopen(12, 15), 3)
    M.put(IV.closedopen(30, 45), 3)
    M.put(IV.closedopen(0, 100), 4)
    print(f"--> Keys: {M.ival}")
    print(f"--> Vals: {M.data}")
    print(f"--> Return: {M.get(IV.closedopen(1, 5))}")

    print("ROUND 2")
    M = IntervalMap()
    M.put(IV.closedopen(30, 45), 1)

    x = IV.closedopen(70, 100)
    y = IV.closedopen(0, 20)
    z = x.union(y)
    M.put(z, 2)

    print(f"--> Keys: {M.ival}")
    print(f"--> Vals: {M.data}")

    M.put(IV.closedopen(15, 35), 3)

    print(f"--> Keys: {M.ival}")
    print(f"--> Vals: {M.data}")

NWK_DIR = "data/synth/nwk"
MTX_DIR = "data/synth/mtx"
PKL_DIR = "data/synth/pkl"
SEQ_DIR = "data/synth/seq"
VIS_DIR = "data/synth/vis"
BLK_DIR = "data/synth/blk"

NITS = 10
MUS  = np.linspace(0, 1e-6, 20)
TS   = np.array([ int(N / n) for n in [7, 5, 3, 1, .5, .25] ], dtype=np.int)

def mk_grid():
    for M, mu in enumerate(MUS):
        for I, T in enumerate(TS):
            for it in range(NITS):
                fa   = f"{SEQ_DIR}/mu{M:02d}_T{I:02d}_{it:02d}.fa"
                bfa  = f"{BLK_DIR}/mu{M:02d}_T{I:02d}_{it:02d}.fa"
                mtx  = f"{MTX_DIR}/mu{M:02d}_T{I:02d}_{it:02d}.mtx"
                tree = f"{NWK_DIR}/mu{M:02d}_T{I:02d}_{it:02d}.nwk"

                print(f"Building random seqs")
                seqs = [random_seq(L, barcode=n) for n in range(N)]
                ancs = seqs
                for _ in range(T):
                    seqs = evolve(seqs, mu)

                # W    = ancestral_weights(seqs)
                print(f"Collecting blocks")
                blks = ancestral_blocks(seqs)

                # Get distribution of block lengths
                lens, nums, blkseqs = [], [], []
                for n, blk in enumerate(blks):
                    for i, iv in enumerate(blk.ival):
                        lens.append(iv.upper-iv.lower)
                        nums.append(len(blk.data[i]))
                        for atom in iv:
                            if atom.upper <= L0:
                                blkseqs.append([ancs[n][0][atom.lower:atom.upper]])
                            else:
                                blkseqs.append([np.array(ancs[n][0][atom.lower:L0].tolist() + \
                                                         ancs[n][0][0:atom.upper].tolist() )])
                lens = np.array(lens)

                print(f"Number of blocks: {len(lens)}")
                # print(f"--> Output {bfa}")
                with open(bfa, "w+") as out:
                    tofasta(blkseqs, out, prefix="block")

                # xratio = np.sum(lens/np.sum(lens) * nums)
                # print(f"Compression Ratio: {xratio}")
                # cdfplot(lens)

                # Output individual current sequences
                # print(f"--> Output {fa}")
                with open(fa, "w+") as out:
                    tofasta(seqs, out)

                # Output kmer guide tree
                mash(fa, mtx)
                # print(f"--> Output {mtx}")
                with open(mtx) as out:
                    D, names = getmtx(out)

                sf  = squareform(D)
                Z   = linkage(sf)
                t   = to_tree(Z)
                nwk = getnwk(t, "", t.dist, names)

                # print(f"--> Output {tree}")
                with open(tree, 'w+') as out:
                    out.write(f"{nwk}\n")

    with open("data/synth/params.npz", "w+") as out:
        np.savez(out, Mus=MUS, Ts=TS)

    # # Run graph creation code
    # f = "test.nwk"
    # Ss = fai.Fasta(fa)

    # g, ok = Graph.fromnwk(f, Ss) #, mu=0, beta=0)
    # g.name = "test"
    # g.tofasta()
    # pickle.dump(g.todict(), open(f"{g.name}.pkl", "wb"))

    # os.system(f"minimap2 -t 2 -x asm5 -D -c blocks.fasta data/graph/aln/test.fasta 1>test.paf 2>log")
    # hits = parsepaf("test.paf")

if __name__ == "__main__":
    # --------------------
    # Possible entry points
    # test_interval_map()
    # mk_grid()
    # --------------------

    import re
    import matplotlib.pylab as plt
    from py.utils import cdfplot

    T  = N // 7
    mu = 0

    STOP = False
    for nit in range(20):
        if STOP:
            break
        print(f"Running iteration {nit:02d}")
        seqs = [random_seq(L, barcode=n) for n in range(N)]
        ancs = seqs.copy()
        for _ in range(T):
            seqs = evolve(seqs, mu)

        print(f"Collecting blocks")
        blks = ancestral_blocks(seqs)

        blkseqs = {}
        for n, blk in enumerate(blks):
            for i, iv in enumerate(blk.ival):
                for j, atom in enumerate(iv):
                    assert atom.lower < L0, "bad interval"
                    if atom.upper < L0:
                        blkseqs[(n, i, j)] = "".join(chr(c) for c in ancs[n][0][atom.lower:atom.upper])
                    else:
                        blkseqs[(n, i, j)] = "".join(chr(c) for c in (ancs[n][0][atom.lower:L0].tolist() + \
                                                                      ancs[n][0][0:atom.upper].tolist()))

        # Convert sequences to strings
        for i, seq in enumerate(seqs):
            seqs[i] = ("".join(chr(c) for c in seq[0]), seq[1], seq[2])

        for (b, i, j), blk in blkseqs.items():
            if len(blk) < 10:
                continue

            locs = []
            for n, seq in enumerate(seqs):
                locs.extend([(n, m.start()) for m in re.finditer(blk, seq[0])])

            data = set()
            for elt in blks[b].data[i]:
                if elt[1][0] < elt[1][1]:
                    if abs(elt[1][1]-elt[1][0] - (len(blk)-1)) < 10:
                        data.add((elt[0], elt[1][0]))
                else:
                    if abs((len(seqs[elt[0]][0]) - elt[1][0]) + elt[1][1] - (len(blk)-1)) < 10:
                        data.add((elt[0], elt[1][1]))

            locs = set(locs)
            if locs != data and len(blk) > 100:
                print(locs-data)
                STOP = True
                break
