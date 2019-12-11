import os
import sys
import io
import traceback
from glob import glob
import json

import numpy as np
# NOTE: Remove this when done debugging...
np.random.seed(42)

from scipy.cluster.hierarchy import dendrogram, linkage, to_tree
from scipy.spatial.distance  import squareform

import pyfaidx as fai
from Bio           import SeqIO, Phylo
from Bio.Seq       import Seq
from Bio.SeqRecord import SeqRecord

from graph import Graph
from util  import parse_paf, partition_cigar
from cigar import Cigar

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

self_maps = 25
# cluster_id     = int(sys.argv[1]) or 1
# clusters_by_id = {int(c.split('_')[-3]):c for c in glob.glob('data/graph/seq/*fa*')}
# cluster        = clusters_by_id[cluster_id][:-6]
# working_dir    = os.path.basename(cluster)+'_dir'
# if not os.path.isdir(working_dir):
#     os.mkdir(working_dir)


# -----------------------------------------------------------------------------
# Functions

def map_and_merge(graph, fname1, fname2, out):
    try:
        os.system(f"minimap2 -x asm5 -D -c  {fname1} {fname2} 1> {out}.paf 2>log")

        paf = parse_paf(f"{out}.paf")
        paf.sort(key=lambda x:-x['aligned_bases'])

        merged_blks = set()
        if len(paf) == 0:
            return graph, False

        merged = False
        for hit in paf:
            if hit['qry']['name'] in merged_blks \
            or hit['ref']['name'] in merged_blks \
            or hit['ref']['name'] == hit['qry']['name']:
                continue

            # NOTE: Not sure why this is here...
            # if set(graph.blks[hit['qry']['name']].muts.keys()).intersection(graph.blks[hit['ref']['name']].muts.keys()):
            #     continue

            cigar_items = list(Cigar(hit['cigar']).items())
            # if np.sum([x[0] for x in cigar_items if x[1]=='M']) < 0.95*hit['aligned_length']:
            #     print("poor match", hit["cigar"])
            #     continue

            print(f"----> Merging {hit['qry']['name']} with {hit['ref']['name']}")
            graph.merge_hit(hit)
            merged = True
            merged_blks.add(hit['ref']['name'])
            merged_blks.add(hit['qry']['name'])

        graph.prune_empty()
        return graph, merged
    except Exception as exc:
        print(f"Panic: {exc}")
        traceback.print_exc()
        import ipdb
        ipdb.set_trace()

def base(name):
    return ".".join(os.path.basename(name).split(".")[:-1])

def torec(seq):
    return SeqRecord(Seq(str(seq)), id=seq.name, name=seq.name)

def all_blks_contained(g, verbose=False):
    Bs = set(g.blks.keys())
    Ss = set([b[0] for s in g.seqs.keys() for b in g.seqs[s]])
    if verbose:
        print(f"All stored blocks: {Bs}")
        print(f"All sequence blocks: {Ss}")
    return Bs == Ss

def compute_suffix_tree_mtx(G, cls):
    # if G.ST is None:
    #     G.make_suffix_tree()
    #     G.strip_suffix_tree()

    # names = list(G.seqs.keys())
    # D     = np.zeros((len(names), len(names)))

    # for n, name1 in enumerate(names):
    #     for m, name2 in enumerate(names[n:]):
    #         D[n, m+n] = len(G.find_common_substrings(set([name1, name2]))[0])
    #         D[m+n, n] = D[n, m+n]

    t     = suffix.Tree(G.seqs)
    names = list(G.seqs.keys())
    D     = np.zeros((len(names), len(names)))

    for n, name1 in enumerate(names):
        for m, name2 in enumerate(names[n:]):
            D[n, m+n] = len(t.matches(name1, name2))
            D[m+n, n] = D[n, m+n]
    json.dump({'mtx': D.tolist(), 'iso': names},
            open(f"{mtxdir}/{cls}.st.json", 'w+'))

def check(seqs, T, G):
    nerror = 0
    uncompressed_length = 0
    for n in T.get_terminals():
        if n.name not in G.seqs:
            continue
        seq = seqs[n.name]
        orig = str(seq[:].seq).upper()
        print(f"Analyzing {n.name}")
        rec = G.extract(n.name)
        uncompressed_length += len(orig)
        print(f"--> Verifying")
        if orig != rec:
            nerror += 1
            # import ipdb; ipdb.set_trace()
            for i in range(len(orig)//100):
                if (orig[i*100:(i+1)*100] != rec[i*100:(i+1)*100]):
                    print("-----------------")
                    print("O:", i, orig[i*100:(i+1)*100])
                    print("G:", i, rec[i*100:(i+1)*100])

                    diffs = [i for i in range(len(rec)) if rec[i] != orig[i]]
                    pos   = [0]
                    blks  = G.seqs[n.name]
                    for b, strand in blks:
                        pos.append(pos[-1] + sum(1 for n in G.blks[b].seq if not n == "-"))
                    pos = pos[1:]

                    import ipdb
                    ipdb.set_trace()
        else:
            print("----> PASS")

    if nerror == 0:
            print("all sequences correctly reconstructed")
            tlength = np.sum([len(x) for x in G.blks.values()])
            print(f"total graph length: {tlength}")
            print(f"total input sequence: {uncompressed_length}")
            print(f"compression: {uncompressed_length/tlength:1.2f}")
    else:
        raise ValueError("bad sequence reconstruction")



def main():
    nwks = glob(f"{nwkdir}/*.nwk")
    seqs = fai.Fasta(allseq)
    for nwk in nwks:
        # if "005" not in nwk:
        #     continue

        print(f"#### PROCESSING {nwk}")

        cls = base(nwk)
        T   = Phylo.read(nwk,'newick')
        if T.count_terminals() == 1:
            continue
        print("kmer tree, total length:", T.total_branch_length())

        print("initializing terminal node graphs")
        for n in T.get_terminals():
            seq     = torec(seqs[n.name])
            n.graph = Graph.from_seq(seq.id, str(seq.seq).upper())
            n.name  = seq.id
            n.fasta_fname = f"{blddir}/{n.name}.fasta"
            n.graph.to_fasta(n.fasta_fname)

        node_count = 0
        print("merging nodes")
        # iso ='GCA_002587885_1_4'
        for n in T.get_nonterminals(order='postorder'):
            node_count += 1
            n.name  = f'NODE_{node_count:07d}'
            print(f" -- node {n.name} with {n.count_terminals()} children")
            n.graph = Graph.fuse(n.clades[0].graph, n.clades[1].graph)

            print("NEW MAPPING ROUND")
            n.graph, _    = map_and_merge(n.graph, n.clades[0].fasta_fname, n.clades[1].fasta_fname, f"{blddir}/{n.name}")
            n.fasta_fname = os.path.join(*[blddir, n.name+'.fasta'])

            contin = True
            i = 0
            while contin:
                n.graph.to_fasta(n.fasta_fname+f'_iter_{i}')
                print(f"--> Merge cycle {i}")
                n.graph, contin = map_and_merge(n.graph, n.fasta_fname+f'_iter_{i}', n.fasta_fname+f'_iter_{i}', f"{blddir}/{n.name}_iter_{i}")
                i += 1

                contin &= i < self_maps

            print(f"  --- Blocks: {len(n.graph.blks)}, length: {np.sum([len(b) for b in n.graph.blks.values()])}")
            n.graph.to_fasta(n.fasta_fname)

            assert all_blks_contained(n.graph), "missing blocks: intermediate"

        G = T.root.graph
        # try:
        #     G.prune_transitive_edges()
        # except:
        #     print("error: could not prune transitive edges")

        assert all_blks_contained(G), "missing blocks"
        G.to_json(f"{visdir}/{cls}.json", min_length=500)
        G.to_fasta(f"{alndir}/{cls}.fa")
        # compute_suffix_tree_mtx(G, cls)

        # check(seqs, T, G)

    return G


# -----------------------------------------------------------------------------
# Main point of entry

if __name__ == "__main__":
    G = main()
