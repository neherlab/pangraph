import os
import argparse, gzip
from enum import IntEnum
from glob import glob
from collections import defaultdict, OrderedDict

import numpy as np
import matplotlib.pylab as plt
import seaborn as sns

import pyfaidx as fai

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

nselfmap = 2

# ------------------------------------------------------------------------
# Functions

def openany(fname, mode = 'r'):
    if fname.endswith('.gz'):
        return gzip.open(fname, mode)
    else:
        return open(fname, mode)

def map_and_merge(graph, fname1, fname2, out):
    os.system(f"minimap2 -x asm5 -D -c  {fname1} {fname2} 1>{out}.paf 2>log")
    paf = parse_paf(f'{out}.paf')
    paf.sort(key=lambda x:-x['aligned_bases'])

    merged_blocks = set()
    for hit in paf:
        if hit['query']['name'] in merged_blocks \
          or hit['ref']['name'] in merged_blocks \
          or hit['ref']['name']==hit['query']['name']:
            continue

        if set(graph.blocks[hit['query']['name']].sequences.keys()).intersection(graph.blocks[hit['ref']['name']].sequences.keys()):
            continue

        cigar_items = list(Cigar(hit['cigar']).items())
        if np.sum([x[0] for x in cigar_items if x[1]=='M']) < 0.95*hit['aligned_length']:
            print("poor match", hit["cigar"])
            continue

        graph.merge_hit(hit)
        merged_blocks.add(hit['ref']['name'])
        merged_blocks.add(hit['query']['name'])
    graph.prune_empty()

    return graph

def partition_pair(qname, Lq, rname, Lr, out):
    os.system(f"minimap2 -x asm5 -D -c  {qname} {rname} 1>{out}.paf 2>log")
    paf = parse_paf(f'{out}.paf')
    assert len(paf) > 0

    paf.sort(key=lambda x:-x['aligned_bases'])
    qry, ref = Partition(Lq), Partition(Lr)
    for hit in paf:
        qry.add_interval(hit['query']['start'], hit['query']['end'])
        ref.add_interval(hit['ref']['start'], hit['ref']['end'])

    return len(qry) + len(ref)

def get_pairs(dirpath):
    hits = defaultdict(lambda: defaultdict(bool))

    for aln in glob(f"{dirpath}/*.paf"):
        ref = os.path.basename(aln)[:-4]
        for line in open(aln, 'r'):
            entry = line.strip().split()
            assert entry[PAF.rname] == ref
            hits[ref][entry[PAF.qname]] = True

    return hits

def main(args):
    # Fasta file must be indexed for random access
    isolates = fai.Fasta(args.fa)

    names = {}
    n     = 0
    for name in isolates.keys():
        if name not in names:
            names[name] = n
            n += 1

    # Compile all hits
    hits = get_pairs(args.dir)
    # with openany(args.dir) as fh:
    #     for line in fh:
    #         entry = line.strip().split()

    #         q, r = str(entry[PAF.qname]).lstrip(">"), str(entry[PAF.rname]).lstrip(">")

    #         # hits[q][r].append( ( (int(entry[PAF.qbeg]), int(entry[PAF.qend]), int(entry[PAF.qlen])), \
    #         #                      (int(entry[PAF.rbeg]), int(entry[PAF.rend]), int(entry[PAF.rlen])), \
    #         #                       int(entry[PAF.alen]), int(entry[PAF.mapq]), entry[PAF.strand] == "+") )

    D = np.inf*np.ones((len(names), len(names)))
    np.fill_diagonal(D, 0)

    # Make graph for each pair. Count number of entries.
    for ref in hits.keys():
        rseq = isolates[ref]
        with open(reffa, 'w') as fh:
            fh.write(f">{str(rseq.name)}\n")
            fh.write(str(rseq[:]))

        ri = names[ref]
        # qG = Graph.from_sequence(qn, qseq[:].seq)
        for qry in hits[ref].keys():
            if ref == qry: #or len(matches) == 0:
                continue

            qi = names[qry]
            if D[qi, ri] != 0:
                D[ri, qi] = D[qi, ri]

            qseq = isolates[qry]
            with open(qryfa, 'w') as fh:
                fh.write(f">{str(qseq.name)}\n")
                fh.write(str(qseq[:]))

            # rG = Graph.from_sequence(rn, rseq[:].seq)
            # fG = Graph.fuse(qG, rG)
            # fG = map_and_merge(fG, qryfa, reffa, mapaf)

            # for _ in range(nselfmap):
            #     fG.to_fasta(outfa)
            #     fG = map_and_merge(fG, outfa, outfa, mapaf)
            # fG.prune_transitive_edges() # This line throws an error

            # D[names[qn],names[rn]] = len(fG.blocks)
            # D[names[rn],names[qn]] = len(fG.blocks)

            D[ri, qi] = partition_pair(reffa, len(rseq[:]), qryfa, len(qseq[:]), mapaf)

    # Serialize the result
    ordered_names = np.empty(len(names), dtype=object)
    for name, n in names.items():
        ordered_names[n] = name

    np.savez("data/graphdist.npz", ordered_names, D)

# ------------------------------------------------------------------------
# Main point of entry

parser = argparse.ArgumentParser(description = "", usage = "get alignment sparsity")
parser.add_argument("fa",  type = str, help = "fasta file")
parser.add_argument("paf", type = str, help = "paf file")

if __name__ == "__main__":
    args = parser.parse_args()
    main(args)
