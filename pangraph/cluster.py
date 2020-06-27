"""
builds a guide tree utilized by the pan-genome alignment
"""
import io
import json
import subprocess as spawn

import numpy as np
import matplotlib.pylab as plt

from Bio import SeqIO

from .tree import Tree

def register_args(parser):
    parser.add_argument("-d", "--dir",
                        metavar="directory",
                        type=str,
                        default=".",
                        help="directory for output file")
    parser.add_argument("-b", "--backend",
                        metavar="backend",
                        type=str,
                        nargs='?',
                        default="mash",
                        choices=['mash'],
                        help="backend used to estimate inter-sequence distance")
    parser.add_argument("input",
                        type=str,
                        nargs='?',
                        default="-",
                        help="fasta file to cluster")

# ------------------------------------------------------------------------
# mash backend

# NOTE: mash takes '-' as filename if it is to read from stdin
def run_mash(inpath):
    stdout = spawn.check_output(f"mash triangle {inpath} 2>/dev/null", shell=True)
    return io.StringIO(stdout.decode("utf-8"))

def parse_mash(input):
    nrows = int(input.readline().strip())
    M     = np.zeros((nrows, nrows), dtype=float)
    names = []
    for i, line in enumerate(input):
        e = line.strip().split()
        names.append(e[0])
        M[i,:i] = [float(x) for x in e[1:]]

    M  = M + M.T

    input.close()
    return M, np.array(names)

# ------------------------------------------------------------------------
# all backends here 

class Backend(object):
    def __init__(self, run, parse):
        self.run   = run
        self.parse = parse

backends = {
    "mash" : Backend(run=run_mash, parse=parse_mash),
    # add more backends here...
}

# gpath = "data/mtx/graphdist.nomap.npz"
# fname = "data/graph/clusters.tsv"
# tdir  = "data/graph/nwk"
# sdir  = "data/graph/seq"
# mdir  = "data/graph/mtx"

# ------------------------------------------------------------------------
# helpers

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

def main(args):
    '''
    Parameters
    ----------
    args : namespace
        arguments passed in via the command-line from pangraph
    Returns
    -------
    int
        returns 0 for success, 1 for general error
    '''

    backend = backends[args.backend]
    outdir  = args.dir.rstrip("/")

    dist, names = backend.parse(backend.run(args.input))
    tree = Tree.nj(dist, names)

    with open(args.input, 'r') as fd:
        seqs = {s.id : s.seq for s in SeqIO.parse(fd, "fasta")}
    tree.attach(seqs)

    # exports:
    np.savez(f"{outdir}/distmtx.npz", dist=dist, names=names)

    # for now we output both forms of the tree
    with open(f"{outdir}/guide.nwk", "w") as fd:
        tree.write_nwk(fd)

    with open(f"{outdir}/guide.json", "w") as fd:
        tree.write_json(fd)

    return 0

# TODO: Remove assumption of having kpath/gpath already computed.
# if __name__ == "__main__":
#     d = np.load(gpath)

#     # Load in kmer distance and graph distance
#     # Permute rows/columns in kmer to graph
#     Dk, knames = parsekmer(kpath)
#     Dg, gnames = d['arr_1'], d['arr_0']

#     Dg = .5 * (Dg + Dg.T)

#     maps = np.zeros(len(gnames))
#     for i, name in enumerate(gnames):
#         j = np.where(knames == name)[0][0]
#         maps[i] = int(j)

#     Dkp = np.zeros_like(Dg)
#     for i in range(len(maps)):
#         for j in range(len(maps)):
#             Dkp[i,j] = Dk[maps[i], maps[j]]

#     Dk = Dkp

#     # Cluster matrix and export results
#     Dg[Dk < .35]     = np.inf
#     Dg[np.isinf(Dg)] = 500

#     Dsq = ssd.squareform(Dg)
#     Z   = sch.linkage(Dsq, method="average")
#     cls = sch.fcluster(Z, 200, criterion="distance")
#     cls = np.array(cls)

#     # Export final clusters as tsv
#     export(cls, gnames, fname)

#     # Export final clusters as newicks
#     C = np.max(cls)
#     T = Tree(asnwkstr(sch.to_tree(Z), "", 0, gnames))
#     for c in range(C):
#         Tp = T.copy()
#         Tp.prune(gnames[cls == (c+1)])
#         with open(f"{tdir}/cluster_{c:03d}.nwk", "w") as out:
#             out.write((Tp.write()))

#     # Export final clusters as fastas
#     seqs = fai.Fasta("data/seq/all.fasta")
#     for c in range(C):
#         with open(f"{sdir}/cluster_{c:03d}.fa", "w") as out:
#             for name in gnames[cls == (c+1)]:
#                 out.write(f">{name}\n")
#                 out.write(f"{str(seqs[name])}\n")

#     # Export submatrices 
#     for c in range(C):
#         indx = cls == (c+1)
#         Dsub = Dg[indx,:]
#         Dsub = Dsub[:, indx]
#         np.savez(f"{mdir}/cluster_{c:03d}.npz", Dsub, gnames[indx])
#         json.dump({'mtx' : Dsub.tolist(), "iso" : gnames[indx].tolist()}, \
#                 open(f"{mdir}/cluster_{c:03d}.json", "w+"))
