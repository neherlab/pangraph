import os
import argparse, gzip
import pickle as pkl
from enum import IntEnum
from glob import glob
from collections import defaultdict, OrderedDict

import numpy as np
import matplotlib.pylab as plt
import seaborn as sns

from Bio import SeqIO

from graph import Graph
from cigar import Cigar
from util  import parse_paf
from interval import Partition

# ------------------------------------------------------------------------
# Global constants/variables

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

qryfa = "tmp/qry.fa"
reffa = "tmp/ref.fa"
outfa = "tmp/out.fa"
mapaf = "tmp/aln"
remap = False

# ------------------------------------------------------------------------
# Functions

def openany(fname, mode = 'r'):
    if fname.endswith('.gz'):
        return gzip.open(fname, mode)
    else:
        return open(fname, mode)

def subselect(paf, qryname):
    return [hit for hit in paf if hit['query']['name'] == qryname]

def partition_pair_w_map(qname, Lq, rname, Lr, out):
    os.system(f"minimap2 -x asm5 -D -c  {qname} {rname} 1>{out}.paf 2>log")
    paf = parse_paf(f'{out}.paf')
    assert len(paf) > 0

    paf.sort(key = lambda x: -x['aligned_bases'])
    qry, ref = Partition(Lq), Partition(Lr)
    for hit in paf:
        qry.add_interval(hit['query']['start'], hit['query']['end'])
        ref.add_interval(hit['ref']['start'], hit['ref']['end'])

    return len(qry) + len(ref), (qry, ref)

def partition_pair(qname, Lq, rname, Lr, paf):
    assert len(paf) > 0

    paf.sort(key = lambda x: -x['aligned_bases'])
    qry, ref = Partition(Lq), Partition(Lr)
    for hit in paf:
        qry.add_interval(hit['query']['start'], hit['query']['end'])
        ref.add_interval(hit['ref']['start'], hit['ref']['end'])

    return len(qry) + len(ref), (qry, ref)

def get_pairs(dirpath):
    hits = defaultdict(lambda: set())

    for aln in glob(f"{dirpath}/*.paf"):
        ref = os.path.basename(aln)[:-4]
        for line in open(aln, 'r'):
            entry = line.strip().split()
            assert entry[PAF.rname] == ref
            hits[ref].add(entry[PAF.qname])

    return hits

def main(args):
    isolates = []
    seqlen   = {}
    for path in glob(f"{args.seqdir}/*.fasta"):
        isolates.append(os.path.basename(path)[:-6])
        seqlen[isolates[-1]] = len(next(SeqIO.parse(path, 'fasta')))

    names = {}
    n     = 0
    for name in sorted(isolates):
        if name not in names:
            names[name] = n
            n += 1

    # Compile all hits
    print("Getting pairs")
    hits = get_pairs(args.alndir)
    print("Done getting pairs")

    D  = np.inf*np.ones((len(names), len(names)))
    np.fill_diagonal(D, 0)
    DB = defaultdict(lambda: defaultdict(tuple))

    for ref in hits.keys():
        ri = names[ref]
        print(f"At isolate {ref}")
        if not remap:
            paf = parse_paf(f"{args.alndir}/{ref}.paf")

        for qry in hits[ref]:
            if ref == qry or qry not in hits: # Second term is temporary
                continue

            qi = names[qry]
            if D[qi, ri] != 0:
                D[ri, qi] = D[qi, ri]

            if remap:
                D[ri, qi], DB[ri][qi] = partition_pair_w_map(reffa, seqlen[ref], qryfa, seqlen[qry], mapaf)
            else:
                D[ri, qi], DB[ri][qi] = partition_pair(reffa, seqlen[ref], qryfa, seqlen[qry], subselect(paf, qry))

    # Serialize the result
    ordered_names = np.empty(len(names), dtype=object)
    for name, n in names.items():
        ordered_names[n] = name

    for key in DB.keys():
        DB[key] = dict(DB[key])
    DB = dict(DB)

    np.savez("data/graphdist.noremap.npz", ordered_names, D, DB)
    pkl.dump(DB, "data/partitions.pkl")

# ------------------------------------------------------------------------
# Main point of entry

parser = argparse.ArgumentParser(description = "", usage = "get alignment sparsity")
parser.add_argument("seqdir", type = str, help = "directory to singleton sequences")
parser.add_argument("alndir", type = str, help = "directory to aligned singletons")

if __name__ == "__main__":
    args = parser.parse_args()
    main(args)
