import argparse, gzip
from enum import IntEnum
from collections import defaultdict

import numpy as np
import scipy.sparse as ssp
import matplotlib.pylab as plt
import seaborn as sns

class PAF(IntEnum):
    qname  = 0
    qlen   = 1
    qbeg   = 2
    qend   = 3
    strand = 4
    rname  = 5
    rlen   = 6
    rbeg   = 7
    rend   = 8
    nmatch = 9
    alen   = 10
    mapq   = 11

def openany(fname, mode = 'r'):
    if fname.endswith('.gz'):
        return gzip.open(fname, mode)
    else:
        return open(fname, mode)

def main(args):
    names = {}

    # Build sparse matrix
    i = 0
    with openany(args.fa, 'rt') as fh:
        for line in fh:
            if line.startswith('>'):
                name = str(line[1:].strip().split()[0])
                if name not in names:
                    names[name] = i
                    i += 1

    hits = ssp.dok_matrix((len(names),len(names)))

    with openany(args.paf) as fh:
        for line in fh:
            entry = line.strip().split()
            q, r = str(entry[PAF.qname]).lstrip(">"), str(entry[PAF.rname]).lstrip(">")
            qname, rname = names[q], names[r]
            hits[qname, rname] += 1

    # Plot sparse matrix
    hits = hits.todense()
    hits[hits==0] = np.nan

    plt.matshow(hits)
    plt.clim(0, 5)
    plt.colorbar()
    plt.show()

    return hits
    # n = np.ndarray.flatten(np.sum(hits > 0, axis=0)[:]) / 2
    # plt.plot(sorted(n), np.linspace(0, 1, len(n)))
    # plt.show()

# ------------------------------------------------------------------------
# Main point of entry

parser = argparse.ArgumentParser(description = "", usage = "get alignment sparsity")
parser.add_argument("paf", type = str, help = "paf file")
parser.add_argument("fa",  type = str, help = "fasta file")

if __name__ == "__main__":
    args = parser.parse_args()
    hits = main(args)
