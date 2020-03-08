import json
import numpy as np
import scipy.spatial.distance as ssd
import scipy.cluster.hierarchy as sch

import pyfaidx as fai
from   ete3 import Tree

# ------------------------------------------------------------------------
# Global constants/variables

kpath = "data/kmerdist.txt"
gpath = "data/mtx/graphdist.nomap.npz"
fname = "data/graph/clusters.tsv"
tdir  = "data/graph/nwk"
sdir  = "data/graph/seq"
mdir  = "data/graph/mtx"

# ------------------------------------------------------------------------
# Helper functions

def parsekmer(mtx):
    with open(mtx) as fh:
        nrows = int(fh.readline().strip())
        M     = np.zeros((nrows, nrows), dtype=float)
        names = []
        for li, line in enumerate(fh):
                e = line.strip().split()
                names.append(e[0].split('/')[-1][:-3])
                M[li,:(li+1)] = [float(x) for x in e[1:]]

        Msym = M + M.T
        Msym -= np.eye(nrows)

    return Msym, np.array(names)

def asnwkstr(node, nwk, pdist, leafs):
    if node.is_leaf():
        return "%s:%.2f%s" % (leafs[node.id], pdist - node.dist, nwk)

def export(cls, names, fname):
    C = np.max(cls)
    with open(fname, "w+") as out:
        for c in range(C):
            out.write("\t".join(names[cls == (c+1)]) + "\n")

# ------------------------------------------------------------------------
# Main point of entry

# TODO: Remove assumption of having kpath/gpath already computed.
if __name__ == "__main__":
    d = np.load(gpath)

    # Load in kmer distance and graph distance
    # Permute rows/columns in kmer to graph
    Dk, knames = parsekmer(kpath)
    Dg, gnames = d['arr_1'], d['arr_0']

    Dg = .5 * (Dg + Dg.T)

    maps = np.zeros(len(gnames))
    for i, name in enumerate(gnames):
        j = np.where(knames == name)[0][0]
        maps[i] = int(j)

    Dkp = np.zeros_like(Dg)
    for i in range(len(maps)):
        for j in range(len(maps)):
            Dkp[i,j] = Dk[maps[i], maps[j]]

    Dk = Dkp

    # Cluster matrix and export results
    Dg[Dk < .35]    = np.inf
    Dg[np.isinf(Dg)] = 500

    Dsq = ssd.squareform(Dg)
    Z   = sch.linkage(Dsq, method="average")
    cls = sch.fcluster(Z, 200, criterion="distance")
    cls = np.array(cls)

    # Export final clusters as tsv
    export(cls, gnames, fname)

    # Export final clusters as newicks
    C = np.max(cls)
    T = Tree(asnwkstr(sch.to_tree(Z), "", 0, gnames))
    for c in range(C):
        Tp = T.copy()
        Tp.prune(gnames[cls == (c+1)])
        with open(f"{tdir}/cluster_{c:03d}.nwk", "w") as out:
            out.write((Tp.write()))

    # Export final clusters as fastas
    seqs = fai.Fasta("data/seq/all.fasta")
    for c in range(C):
        with open(f"{sdir}/cluster_{c:03d}.fa", "w") as out:
            for name in gnames[cls == (c+1)]:
                out.write(f">{name}\n")
                out.write(f"{str(seqs[name])}\n")

    # Export submatrices 
    for c in range(C):
        indx = cls == (c+1)
        Dsub = Dg[indx,:]
        Dsub = Dsub[:, indx]
        np.savez(f"{mdir}/cluster_{c:03d}.npz", Dsub, gnames[indx])
        json.dump({'mtx' : Dsub.tolist(), "iso" : gnames[indx].tolist()}, \
                open(f"{mdir}/cluster_{c:03d}.json", "w+"))
