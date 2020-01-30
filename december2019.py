import simplesam as ssam
import pyfaidx as fai
import numpy as np

import matplotlib.pylab as plt

from glob import glob
from collections import defaultdict

# -----------------------------------------------------------------------
# Global constants

ROOT_1 = "data/graph/aln"
ROOT_2 = "data/graph_2/aln"

# -----------------------------------------------------------------------
# Utility functions

def panic(msg):
    raise ValueError(f"Panic: {msg}")

def getstats(fldr):
    fracmap, alnlen, nummap = [], [], []
    for sam in glob(fldr):
        aln = ssam.Reader(open(sam, "r"))
        seq = fai.Fasta(sam.replace("sam", "fa"))

        print(f"Analyzing {sam}")

        mapfrc = defaultdict(lambda: 0)
        mapnum = defaultdict(lambda: 0)
        for read in aln:
            qry = read.qname
            ref = read.rname
            if read.rname == "*":
                continue

            if qry is not None and ref is not None:
                qseq = seq[qry][:].seq
                rseq = seq[ref][:].seq
                alen = sum(c[0] for c in read.cigars if c[1] == 'M')
                mapfrc[qry] += alen
                mapnum[qry] += 1
                if not qry == ref:
                    mapfrc[ref] += alen
                    mapnum[ref] += 1

        for i, l in mapfrc.items():
            mapfrc[i] = l/len(seq[i][:].seq)

        data = [(mapfrc[i], mapnum[i], len(seq[i][:].seq)) for i in mapfrc.keys()]
        x, y, z = zip(*data)
        fracmap.extend(x)
        nummap.extend(y)
        alnlen.extend(z)

    nummap, alnlen, fracmap = np.array(nummap), np.array(alnlen), np.array(fracmap)

    nummap  = nummap[fracmap > .8]
    alnlen  = alnlen[fracmap > .8]
    fracmap = fracmap[fracmap > .8]

    return np.array(fracmap), np.array(nummap), np.array(alnlen)

if __name__ == "__main__":
    f1, n1, a1 = getstats(f"{ROOT_1}/*.sam")
    f2, n2, a2 = getstats(f"{ROOT_2}/*.sam")

    plt.scatter(f1, a1, s=n1, alpha=1)
    plt.scatter(f2, a2, s=n2, alpha=.25)

    plt.yscale("log")
    plt.xlabel("Fractional alignment degeneracy")
    plt.ylabel("Block length (bp)")
    plt.title("Size = number of hits")
