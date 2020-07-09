"""
build a pangenome alignment from an annotated guide tree
"""
import os, sys
import builtins

from .utils import mkdir, log
from .tree import Tree

def open(path, *args, **kwargs):
    if path == '-':
        return sys.stdin
    return builtins.open(path, *args, **kwargs)

def register_args(parser):
    parser.add_argument("-d", "--dir",
                        metavar="directory",
                        type=str,
                        default=".",
                        help="directory used for output files")
    parser.add_argument("-l", "--len",
                        metavar="cutoff length",
                        type=int,
                        default=50,
                        help="minimum block size for nucleotides")
    parser.add_argument("-m", "--mu",
                        metavar="block cut chemical potential",
                        type=int,
                        default=100,
                        help="energy cost for cutting blocks (used during block merges)")
    parser.add_argument("-x", "--extensive",
                        default=False,
                        action='store_true',
                        help="boolean flag that toggles whether the energy associated to a cut grows proportionally to the number of sequences cut")
    parser.add_argument("-b", "--beta",
                        metavar="block diversity cost",
                        type=int,
                        default=22,
                        help="energy cost for mutations (used during block merges)")
    parser.add_argument("input",
                        type=str,
                        default="-",
                        help="input guide tree [json]")

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
    with open(args.input, 'r') as input:
        T = Tree.from_json(input)

    root = args.dir.rstrip('/')
    tmp  = f"{root}/tmp"
    i    = 0
    while os.path.isdir(tmp) and i < 32:
        i += 1
        tmp = f"{root}/tmp{i:03d}"
    mkdir(tmp)
    T.align(tmp, args.len, args.mu, args.beta, args.extensive)
    # TODO: when debugging phase is done, remove tmp directory

    graphs = T.collect()

    for i, g in enumerate(graphs):
        log(f"graph {i}: nseqs: {len(g.seqs)} nblks: {len(g.blks)}")

    for i, g in enumerate(graphs):
        with open(f"{root}/graph_{i:03d}.fa", 'w') as fd:
            g.write_fasta(fd)

    T.write_json(sys.stdout, no_seqs=True)
    # with open(f"{root}/pangraph.json", "w") as fd:
    #     T.write_json(fd, no_seqs=True)

    return 0
