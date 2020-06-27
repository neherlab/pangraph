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
# backends perform the pairwise distance approximation

class Backend(object):
    def __init__(self, run, parse):
        self.run   = run
        self.parse = parse

backends = {
    "mash" : Backend(run=run_mash, parse=parse_mash),
    # add more backends here...
}

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
