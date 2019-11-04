import os
import sys
import io
from glob import glob
import json

import numpy as np
from scipy.cluster.hierarchy import dendrogram, linkage, to_tree
from scipy.spatial.distance  import squareform

import pyfaidx as fai
from Bio           import SeqIO, Phylo
from Bio.Seq       import Seq
from Bio.SeqRecord import SeqRecord

from graph     import Graph
from util      import parse_paf
from cigar     import Cigar

import suffix

# -----------------------------------------------------------------------------
# Global constants/parameters

allseq = "data/seq/all.fasta"
seqdir = "data/seq/singles"
outdir = "data/graph"
visdir = f"{outdir}/vis"
alndir = f"{outdir}/aln"
nwkdir = f"{outdir}/nwk"
blddir = f"{outdir}/bld"
mtxdir = f"{outdir}/mtx"

self_maps = 2
# cluster_id     = int(sys.argv[1]) or 1
# clusters_by_id = {int(c.split('_')[-3]):c for c in glob.glob('data/graph/seq/*fa*')}
# cluster        = clusters_by_id[cluster_id][:-6]
# working_dir    = os.path.basename(cluster)+'_dir'
# if not os.path.isdir(working_dir):
#     os.mkdir(working_dir)


# -----------------------------------------------------------------------------
# Functions

def map_and_merge(graph, fname1, fname2, out):
    os.system(f"minimap2 -x asm5 -D -c  {fname1} {fname2} 1> {out}.paf 2>log")

    paf = parse_paf(f"{out}.paf")
    paf.sort(key=lambda x:-x['aligned_bases'])

    merged_blocks = set()
    for hit in paf:
        if hit['qry']['name'] in merged_blocks \
        or hit['ref']['name'] in merged_blocks \
        or hit['ref']['name'] == hit['qry']['name']:
            continue

        if set(graph.blocks[hit['qry']['name']].sequences.keys()).intersection(graph.blocks[hit['ref']['name']].sequences.keys()):
            continue

        cigar_items = list(Cigar(hit['cigar']).items())
        if np.sum([x[0] for x in cigar_items if x[1]=='M']) < 0.95*hit['aligned_length']:
            print("poor match", hit["cigar"])
            continue

        graph.merge_hit(hit)
        merged_blocks.add(hit['ref']['name'])
        merged_blocks.add(hit['qry']['name'])

    graph.prune_empty()
    return graph

def base(name):
    return ".".join(os.path.basename(name).split(".")[:-1])

def torec(seq):
    return SeqRecord(Seq(str(seq)), id=seq.name, name=seq.name)

def compute_suffix_tree_mtx(G, cls):
    # if G.ST is None:
    #     G.make_suffix_tree()
    #     G.strip_suffix_tree()

    # names = list(G.sequences.keys())
    # D     = np.zeros((len(names), len(names)))

    # for n, name1 in enumerate(names):
    #     for m, name2 in enumerate(names[n:]):
    #         D[n, m+n] = len(G.find_common_substrings(set([name1, name2]))[0])
    #         D[m+n, n] = D[n, m+n]

    t = suffix.Tree(G.sequences)
    names = list(G.sequences.keys())
    D     = np.zeros((len(names), len(names)))

    for n, name1 in enumerate(names):
        for m, name2 in enumerate(names[n:]):
            D[n, m+n] = len(t.matches(name1, name2))
            D[m+n, n] = D[n, m+n]
    json.dump({'mtx': D.tolist(), 'iso': names},
            open(f"{mtxdir}/{cls}.st.json", 'w+'))

def main():
    nwks = glob(f"{nwkdir}/*.nwk")
    seqs = fai.Fasta(allseq)

    divs = []
    for nwk in nwks:
        try:
            print(f"processing {nwk}")
            cls = base(nwk)
            T   = Phylo.read(nwk,'newick')
            if T.count_terminals() == 1:
                continue
            print("kmer tree, total length:", T.total_branch_length())

            print("initializing terminal node graphs")
            for n in T.get_terminals():
                seq     = torec(seqs[n.name])
                n.graph = Graph.from_sequence(seq.id, str(seq.seq).upper())
                n.name  = seq.id
                n.fasta_fname = f"{blddir}/{n.name}.fasta"
                n.graph.to_fasta(n.fasta_fname)

            node_count = 0
            print("merging nodes")
            for n in T.get_nonterminals(order='postorder'):
                node_count += 1
                n.name      = f'NODE_{node_count:07d}'
                print(f" -- node {n.name} with {n.count_terminals()} children")

                n.graph = Graph.fuse(n.clades[0].graph, n.clades[1].graph)
                n.graph = map_and_merge(n.graph, n.clades[0].fasta_fname, n.clades[1].fasta_fname, f"{blddir}/{n.name}")
                n.fasta_fname = os.path.join(*[blddir, n.name+'.fasta'])

                for i in range(self_maps):
                    print(f"   --- self map")
                    n.graph.to_fasta(n.fasta_fname+f'_iter_{i}')
                    n.graph = map_and_merge(n.graph, n.fasta_fname+f'_iter_{i}', n.fasta_fname+f'_iter_{i}', f"{blddir}/{n.name}_iter_{i}")
                    print(f"   --- Blocks: {len(n.graph.blocks)}, length: {np.sum([len(b) for b in n.graph.blocks.values()])}")

                print(f"  --- Blocks: {len(n.graph.blocks)}, length: {np.sum([len(b) for b in n.graph.blocks.values()])}")
                n.graph.to_fasta(n.fasta_fname)

            G = T.root.graph
            try:
                G.prune_transitive_edges()
            except:
                print("error: could not prune")

            nnn = 0
            for b in G.blocks.values():
                if nnn > 20:
                    break
                _seqs = b.sequences
                N    = len(_seqs)
                if N < 10:
                    continue
                L   = len(b.consensus)
                # if L < 1000:
                #     continue
                # div = sum( sum(1 for mut in muts.values() if mut != '-') for muts in _seqs.values()) / N
                divs.append(L)
                nnn += 1
        except:
            print("error")

    return divs

# -----------------------------------------------------------------------------
# Main point of entry

if __name__ == "__main__":
    divs = main()
