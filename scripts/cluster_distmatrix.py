import numpy as np
import scipy.spatial.distance as ssd
import scipy.cluster.hierarchy as sch

import matplotlib as mpl
import matplotlib.pylab as plt
import seaborn as sns

import pyfaidx as fai
from ete3 import Tree

from scripts.util import parse_kmer

# ------------------------------------------------------------------------
# Global constants/variables

gmtx  = "data/graphdist.nomap.npz"
kmtx  = "data/kmerdist.txt"
sfas  = "data/seq/all.fasta"
fname = "data/graph/clusters.tsv"
tdir  = "data/graph/nwk"
sdir  = "data/graph/seq"

plot = False

# ------------------------------------------------------------------------
# Functions

def export(cls, names, fname):
    C = np.max(cls)
    with open(fname, "w+") as out:
        for c in range(C):
            out.write("\t".join(names[cls == (c+1)]) + "\n")

def tonwk(node, nwk, pdist, leafs):
    if node.is_leaf():
        return "%s:%.2f%s" % (leafs[node.id], pdist - node.dist, nwk)
    else:
        if len(nwk) > 0:
            nwk = "):%.2f%s" % (pdist - node.dist, nwk)
        else:
            nwk = ");"
        nwk = tonwk(node.get_left(), nwk, node.dist, leafs)
        nwk = tonwk(node.get_right(), ",%s" % (nwk), node.dist, leafs)
        nwk = "(%s" % (nwk)
        return nwk

def main(D, Dk, names):
    D[Dk < .35]    = np.inf
    D[np.isinf(D)] = 500

    Dsq = ssd.squareform(D)
    Z   = sch.linkage(Dsq, method="average")
    cls = sch.fcluster(Z, 200, criterion="distance")
    cls = np.array(cls)

    # Export final clusters as tsv
    export(cls, names, fname)

    # Export final clusters as newicks
    C = np.max(cls)
    T = Tree(tonwk(sch.to_tree(Z), "", 0, names))
    for c in range(C):
        Tp = T.copy()
        Tp.prune(names[cls == (c+1)])
        with open(f"{tdir}/cluster_{c:03d}.nwk", 'w') as out:
            out.write((Tp.write()))

    # Export final clusters as fastas
    seqs = fai.Fasta(sfas)
    for c in range(C):
        with open(f"{sdir}/cluster_{c:03d}.fa", 'w') as out:
            for name in names[cls == (c+1)]:
                out.write(f">{name}\n")
                out.write(f"{str(seqs[name])}\n")

    # Debugging
    if plot:
        C      = np.max(cls)
        colors = mpl.cm.jet(np.linspace(0, 1, C))
        colors = np.array([colors[c, :] for c in cls])
        sns.clustermap(D, method="average", vmin=0, vmax=100, row_colors=colors)
        plt.show()

if __name__ == "__main__":
    # Preprocess arguments
    data = np.load(gmtx)
    Dgraph, graph_names = data['arr_1'], data['arr_0']
    Dkmers, kmers_names = parse_kmer(kmtx)

    Dgraph = .5 * (Dgraph + Dgraph.T)

    mapping = np.zeros(len(graph_names))
    for n, name in enumerate(graph_names):
        i = np.where(kmers_names == name)[0][0]
        mapping[n] = int(i)

    Dkmer_t = np.zeros(Dgraph.shape)
    for i in range(len(mapping)):
        for j in range(len(mapping)):
            Dkmer_t[i, j] = Dkmers[int(mapping[i]), int(mapping[j])]

    Dkmers = Dkmer_t

    # Run code
    main(Dgraph, Dkmers, graph_names)

