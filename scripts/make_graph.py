import os
import sys
import io
from glob import glob

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

# -----------------------------------------------------------------------------
# Global constants/parameters

allseq = "data/seq/all.fasta"
seqdir = "data/seq/singles"
outdir = "data/graph"
nwkdir = f"{outdir}/nwk"
blddir = f"{outdir}/bld"

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

def main():
    nwks = glob(f"{nwkdir}/*.nwk")
    seqs = fai.Fasta(allseq)
    for nwk in nwks:
        # if "101" not in nwk:
        #     continue

        print(f"processing {nwk}")
        cls = base(nwk)
        T   = Phylo.read(nwk,'newick')
        if T.count_terminals() == 1:
            continue
        print("kmer tree, total length:", T.total_branch_length())
        # Phylo.draw_ascii(T)

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
            raise ValueError("panic")
        # To fasta here
        G.to_json(f"{outdir}/{cls}.json")

    return G

# -----------------------------------------------------------------------------
# Main point of entry

if __name__ == "__main__":
    G = main()

# ni = 0
# for seq in seqs:
# 	seq.id += f"_{ni:03d}"
# 	ni+=1

# nerror = 0
# uncompressed_length = 0
# G = T.root.graph
# G.prune_transitive_edges()

# G.make_suffix_tree()
# seqs = ['GCA_002082705_1_2', 'GCA_001806485_1_2']
# sequences = set(seqs)
# cs, covered = G.find_common_substrings(sequences)

# from itertools import combinations
# for (s1,s2), (i1,i2) in zip(combinations(G.sequences, r=2),combinations(range(len(seqs)), r=2)):
# 	# print("\n", s1, s2,"\n")
# 	cs, covered = G.find_common_substrings(set([s1,s2]))
# 	print(s1,s2,covered)

# check_sequences = False
# export_json = False
# subgraph_matrix = False

# if check_sequences:
#     for n in T.get_terminals():
#             seq = seqs[n.name]
#             orig = str(seq.seq).upper()
#             rec = G.extract(n.name)
#             uncompressed_length += len(orig)
#             if (orig!=rec):
#                     nerror += 1
#                     import ipdb; ipdb.set_trace()
#                     for i in range(len(orig)//100):
#                             if (orig[i*100:(i+1)*100]!=rec[i*100:(i+1)*100]):
#                                     # for j,o,r in zip(range(100), orig[i*100:(i+1)*100], rec[i*100:(i+1)*100]):
#                                     # 	if o!=r:
#                                     # 		print(i*100+j,o,r)
#                                     # node_path = T.get_path(n.name)
#                                     # for ni,pn in enumerate(node_path):
#                                     # 	print(ni, pn.name, orig==pn.graph.extract(n.name))
#                                     # import ipdb; ipdb.set_trace()
#                                     print(i,orig[i*100:(i+1)*100])
#                                     print(i,rec[i*100:(i+1)*100])

#     if nerror==0:
#             print("all sequences correctly reconstructed")
#             tlength = np.sum([len(x) for x in G.blocks.values()])
#             print(f"total graph length: {tlength}")
#             print(f"total input sequence: {uncompressed_length}")
#             print(f"compression: {uncompressed_length/tlength:1.2f}")

# if export_json:
# 	T.root.graph.to_json(os.path.join(blddir, f'graph_{cluster_id:03d}.json'))

## make subgraphs
#
# if subgraph_matrix:
#     from itertools import combinations
#     total_errors = 0
#     DM = np.zeros((len(seqs), len(seqs)))
#     for (s1,s2), (i1,i2) in zip(combinations(G.sequences, r=2),combinations(range(len(seqs)), r=2)):
#             # print("\n", s1, s2,"\n")
#             S = G.sub_graph((s1,s2))
#             S.prune_transitive_edges()
#             for s in S.sequences:
#                     seq = seqs[s]
#                     orig = str(seq.seq).upper()
#                     rec = S.extract(s)
#                     if (orig!=rec):
#                             nerror += 1
#             DM[i1,i2] = len(S.blocks)
#             DM[i2,i1] = len(S.blocks)
#             total_errors += nerror
#             if nerror:
#                     print("reconstruction error for pair",(s1,s2), nerror)
#                     import ipdb;ipdb.set_trace()

#     if total_errors:
#             print("there were reconstruction errors")
#     else:
#             print("all pairwise subgraphs led to correct reconstructions.")


#     condensed_DM = squareform(DM)
#     Z = linkage(condensed_DM)
#     Tgraph = to_tree(Z)
#     nwk = getNewick(Tgraph, "", Tgraph.dist, list(G.sequences))
#     Tgraph2 = Phylo.read(io.StringIO(nwk),'newick')
#     Tgraph2.ladderize()

#     Phylo.draw_ascii(Tgraph2)
