import pickle
import gzip

from copy import deepcopy
from collections import defaultdict

import portion as P

import numpy as np
import numpy.random as rng
# rng.seed(0)

# ------------------------------------------------------------------------
# utility functions

def transform(interval, domain, func):
    fwd = lambda I: P.closedopen(I.lower-domain.lower+interval.lower, I.upper-domain.upper+interval.upper)
    rev = lambda I: P.closedopen(I.lower-interval.lower+domain.lower, I.upper-interval.upper+domain.upper)

    return P.Interval(*[fwd(new_atom) for atom in interval for new_atom in func(rev(atom)) if not new_atom.empty])

def bounds(interval):
    return I.upper - I.lower

def domain(intervals):
    I = P.empty()
    for key in intervals.keys():
        I = I.union(key)
    return I

def area(interval):
    return sum(atom.upper - atom.lower for atom in interval)

# ------------------------------------------------------------------------
# main data structures

class Ancestor:
    def __init__(self, anc, iv):
        self.a = anc
        self.i = iv

    def __str__(self):
        return f"{self.a}-{self.i}"

    def __repr__(self):
        return str(self)

    @property
    def area(self):
        return area(self.i)

    def transform(self, domain, func):
        return Ancestor(self.a, transform(self.i, domain, func))

class Sequence:
    def __init__(self, intervals=[], ancestors=[]):
        self.intervals = intervals
        self.ancestors = ancestors

    def __str__(self):
        return str(dict(zip(self.intervals,self.ancestors)))

    def __repr__(self):
        return str(self)

    def __len__(self):
        return len(self.intervals)

    @classmethod
    def ancestor(cls, n, L):
        intervals = [P.closedopen(0, L)]
        ancestors = [Ancestor(n, P.closedopen(0,L))]
        return cls(intervals, ancestors)

    @classmethod
    def splice(cls, S1, I1, S2, I2):
        intervals = []
        ancestors = []

        fwd = lambda I: P.closedopen(I.lower-I1.lower+I2.lower, I.upper-I1.upper+I2.upper)
        rev = lambda I: P.closedopen(I.lower-I2.lower+I1.lower, I.upper-I2.upper+I1.upper)

        for I, a in S1.items():
            if I.overlaps(I1):
                d, o = I-I1, I&I1
                if not d.empty:
                    intervals.append(d)
                    ancestors.append(a.transform(I,lambda I: I-I1))

                # import ipdb; ipdb.set_trace()
                for II, aa in S2[P.Interval(*[fwd(atom) for atom in o])].items():
                    intervals.append(P.Interval(*[rev(atom) for atom in II]))
                    ancestors.append(aa)
            else:
                intervals.append(I)
                ancestors.append(a)

        return cls(intervals, ancestors)

    @classmethod
    def copy(cls, seq):
        return cls(deepcopy(seq.intervals), deepcopy(seq.ancestors))

    @property
    def domain(self):
        I = P.empty()
        for ia in self.intervals:
            I = I.union(ia)
        return I

    def items(self):
        for i, v in zip(self.intervals, self.ancestors):
            yield i, v

    def keys(self):
        for i in self.intervals:
            yield i

    def values(self):
        for v in self.ancestors:
            yield v

    def __getitem__(self, key):
        if isinstance(key, P.Interval):
            keys = []
            vals = []
            for i, v in self.items():
                if i.overlaps(key):
                    keys.append(i & key)
                    vals.append(v.transform(i, lambda I: I & key))
            return Sequence(keys, vals)
        else:
            raise ValueError(f"type {type(key)} not handled yet")

    def __setitem__(self, key, val):
        if isinstance(key, P.Interval):
            rm_index  = []
            new_keys  = []
            new_vals  = []
            for i, (I, v) in enumerate(self.items()):
                if I.overlaps(key):
                    remain = I - key
                    tmp    = key - I

                    rm_index.append(i)
                    if not remain.empty:
                        new_keys.append(remain)
                        new_vals.append(v.transform(I, lambda I: I-key))

                    val, tmp = val.transform(key, lambda I: key-I), key

            if not key.empty:
                new_keys.append(key)
                new_vals.append(val)

            for i in sorted(rm_index, reverse=True):
                del self.intervals[i]
                del self.ancestors[i]

            self.intervals.extend(new_keys)
            self.ancestors.extend(new_vals)
        else:
            raise ValueError(f"type {type(key)} not handled yet")

# ------------------------------------------------------------------------
# simulation

def horizontal_transfer(S1, S2):
    D1, D2 = S1.domain, S2.domain

    lb1 = rng.randint(D1.lower, high=D1.upper)
    ub1 = rng.randint(lb1, high=D1.upper)
    I1  = P.closedopen(lb1, ub1)

    delta = ub1 - lb1

    lb2 = rng.randint(D2.lower, (D2.upper-delta))
    ub2 = lb2 + delta
    I2  = P.closedopen(lb2, ub2)

    return I1, I2

def evolver(r):
    def evolve(seqs):
        rand   = rng.rand(len(seqs))
        parent = rng.randint(len(seqs), size=len(seqs))

        def offspring(i):
            if rand[i] < r:
                S1, S2 = seqs[parent[i]], seqs[rng.randint(len(seqs))]
                I1, I2 = horizontal_transfer(S1, S2)
                seq = Sequence.splice(S1, I1, S2, I2)
                return seq
            else:
                return Sequence.copy(seqs[parent[i]])

        return [offspring(i) for i in range(len(seqs))]

    return evolve

# return the "true" partioning
def tiles(seqs):
    tiles = { i: P.IntervalDict() for i in range(len(seqs)) }
    union = lambda A, B: A.union(B)

    # get the atomic partitions
    for i, seq in enumerate(seqs):
        print(f"sequence {i}")
        for j, (I, anc) in enumerate(seq.items()):
            tiles[anc.a] = tiles[anc.a].combine(P.IntervalDict({anc.i:set([(i,j)])}), union)

    empty = []
    for k, v in tiles.items():
        if len(v) == 0:
            empty.append(k)

    for k in empty:
        tiles.pop(k)

    shift  = lambda iv, old, new: P.closedopen(iv.lower-old.lower+new.lower, iv.upper-old.upper+new.upper)
    tiling = defaultdict(set)
    for i, seq in enumerate(seqs):
        for (I, anc) in seq.items():
            for ai in tiles[anc.a][anc.i].keys():
                for atom in ai:
                    tiling[(anc.a,atom)].add((i,I&shift(atom,anc.i,I)))

    for key, tile in tiling.items():
        for member in tile:
            if area(key[1]) != area(member[1]):
                print("ERROR")
                import ipdb; ipdb.set_trace()

    return tiling

# verification
def ok(*seqs):
    for seq in seqs:
        for iv, anc in seq.items():
            if area(iv) != anc.area:
                print(f"--> failure: {iv} ({area(iv)}) <=> {anc.i} ({anc.area})")
                return False

    for seq in seqs:
        if seq.domain != P.closedopen(0,L):
            print(f"--> failure: {seq.domain}")
            return False

    return True

def random_sequence(length):
    nucleotide = np.array(["A", "C", "G", "T"])
    return "".join(nucleotide[i] for i in rng.randint(len(nucleotide), size=length))

def generate(seqs, tiles, L):
    blocks = defaultdict(P.IntervalDict)
    for ancestor, atom in tiles.keys():
        blocks[ancestor][atom] = random_sequence(area(atom))

    for seq in seqs:
        result = np.chararray(L)
        for I in seq.keys():
            for atom in I:
                block = seq[atom]
                assert len(block) == 1, "bad atom"

                tile = next(block.values())
                assert area(atom) == tile.area, "unequal lengths"

                anc = blocks[tile.a][tile.i]
                dom = domain(anc)
                fn  = lambda I: P.closedopen(I.lower - dom.lower + atom.lower, I.upper - dom.upper + atom.upper)
                for i, s in anc.items():
                    ii = fn(i)
                    for l, c in enumerate(s):
                        result[l+ii.lower] = c

        yield result.tostring().decode('utf-8')

def fasta(seqs, io):
    NL = "\n"
    C  = 80
    for i, seq in enumerate(seqs):
        io.write(f">isolate_{i:04d}{NL}")
        io.write(f"{NL.join(seq[(0+i):(C+i)] for i in range(0,len(seq),C))}{NL}")

N = 100
L = 100000
if __name__ == "__main__":
    evolve = evolver(0.025)
    seqs   = [Sequence.ancestor(n, L) for n in range(N)]
    print("simulating...")
    for t in range(100):
        print(f"--> iteration {t}")
        seqs = evolve(seqs)
        if not ok(*seqs):
            print("ERROR")
            import ipdb; ipdb.set_trace()

    print("collecting...")
    tiling = tiles(seqs)

    print("saving...")
    with gzip.open("data/generated/assemblies/isolates.fna.gz", 'wt') as io:
        fasta(generate(seqs, tiling, L), io)

    with open("data/generated/tiling.pkl", 'wb') as io:
        output = {'tiles': tiling, 'seqs': seqs}
        pickle.dump(output, io)
