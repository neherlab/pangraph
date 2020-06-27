"""
generates synthetic data from a recombinative wright-fisher model
"""

import os
import json

from Bio import SeqIO

from .simulation import Population
from .utils import mkdir, panic, asrecord

def register_args(parser):
    parser.add_argument("-d", "--dir",
                        metavar="directory",
                        type=str,
                        default="data/synth",
                        help="directory for output files")
    parser.add_argument("-n", "--name",
                        metavar="name",
                        type=str,
                        default="generated",
                        help="base name for output")
    parser.add_argument("-N", "--size",
                        metavar="size",
                        type=int,
                        default=100,
                        help="simulated population size")
    parser.add_argument("-L", "--len",
                        metavar="length",
                        type=int,
                        default=100000,
                        help="average (and initial) genome length")
    parser.add_argument("-m", "--mu",
                        metavar="mutation rate",
                        type=float,
                        default=1e-5,
                        help="rate of mutation/genome/generation")
    parser.add_argument("--rate_hgt",
                        metavar="hgt rate",
                        type=float,
                        default=0.0,
                        help="rate of horizontal gene transfer/genome/generation")
    parser.add_argument("--rate_indel",
                        metavar="indel rate",
                        type=float,
                        default=0.0,
                        help="rate of indels/genome/generation")
    parser.add_argument("--rate_transpose",
                        metavar="transposition rate",
                        type=float,
                        default=0.0,
                        help="rate of transpositions/genome/generation")
    parser.add_argument("-T", "--gen",
                        metavar="generations",
                        type=int,
                        default=35,
                        help="number of generations to simulate")

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

    args = vars(args)
    time = args.pop('gen')
    out  = args.pop('dir')
    name = args.pop('name')

    # remove unneeded data
    del args['version']
    del args['__command__']

    pop = Population(**args)
    pop.evolve(time)

    if not mkdir(out):
        panic("failed to make output directory")
        return 1

    # outputs:

    # present day sequences
    with open(f"{out}/seq.fa", 'w') as fd:
        pop.write_fasta(fd)

    blk = pop.ancestral_blocks()

    # ancestral blocks
    with open(f"{out}/ancestral.json", 'w') as fd:
        json.dump(blk, fd)

    # ancestral block sequences
    ancs = []
    for a, b in blk.items():
        nm = f"anc_{a:05d}"
        for i, ivs in enumerate(b['interval']):
            for j, iv in enumerate(ivs):
                if (iv[0] < pop.L):
                    seq = "".join(chr(c) for c in pop.anc[a][0][iv[0]:iv[1]])
                else:
                    seq = "".join(chr(c) for c in pop.anc[a][0][iv[1]:pop.L]) + \
                          "".join(chr(c) for c in pop.anc[a][0][0:iv[0]])
            ancs.append(asrecord(seq, f"{nm}_{j:03d}_{i:03d}"))

    with open(f"{out}/ancestral.fa", 'w') as fd:
        SeqIO.write(ancs, fd, "fasta")

    return 0
