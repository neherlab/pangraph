import os, sys, io, json
from glob        import glob
from collections import defaultdict, Iterable

import numpy as np
from scipy.cluster.hierarchy import dendrogram, linkage, to_tree
from scipy.spatial.distance  import squareform

import matplotlib as mpl
import matplotlib.pylab as plt

import pyfaidx as fai
from Bio           import SeqIO, Phylo
from Bio.Seq       import Seq
from Bio.SeqRecord import SeqRecord

import suffix
from graph import Graph
from util  import parse_paf
from cigar import Cigar

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

def clonality_and_length():
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

# ------------------------------------------------------------------------
# Counter dictionary functions

def default_to_normal(d):
    if isinstance(d, defaultdict):
        d = {k : default_to_normal(v) for k, v in d.items()}
    return d

def median(d):
    if not isinstance(d, defaultdict):
        return np.median(d)

    pairs = sorted(list(d.items()))

    csum = 0
    for i, p in enumerate(pairs):
        pairs[i] = (p[0], p[1] + csum)
        csum = p[1]

    pairs = [(p[0], p[1]/csum) for p in pairs]
    val, num = zip(*pairs)
    val = np.array(val)
    num = np.array(num)

    idx = np.argmin(np.abs(num-.5))
    return val[idx]

def merge_dicts(kmer_locs, as_dict=False):
    all_locs = defaultdict(lambda: defaultdict(lambda: 0)) if as_dict else defaultdict(list)
    for db in kmer_locs.values():
        for Z, pairs in db.items():
            if not as_dict:
                all_locs[Z].extend(pairs)
            else:
                for pos, num in pairs.items():
                    all_locs[Z][pos] += num

    return default_to_normal(all_locs)

# ------------------------------------------------------------------------
# Kmer location functions

# Shifting kmers in window
def shiftfn(k):
    mask, cmpl = (1 << (2*k)) - 1, (2 * (k-1))

    def f(c, Z_fwd, Z_rev):
        Z_fwd = ( (Z_fwd << 2) | c ) & mask
        Z_rev = ( (Z_rev >> 2) | ((3 ^ c) << cmpl) )
        return Z_fwd, Z_rev

    return f

# Mapping from nucleotide to integer
nuctab = 4 * np.ones(256, dtype=np.int)

nuctab[ord('a')] = 0
nuctab[ord('A')] = 0
nuctab[ord('c')] = 1
nuctab[ord('C')] = 1
nuctab[ord('g')] = 2
nuctab[ord('G')] = 2
nuctab[ord('t')] = 3
nuctab[ord('T')] = 3

def bin_kmer_location(k=8, normalize=True, consensus=True, as_dict=False):
    kmer_locs = {}

    shift = shiftfn(k)
    n = 0
    for path in glob(f"{visdir}/*.json"):
        print(f"----> Analyzing {path}")
        locs = defaultdict(lambda: defaultdict(lambda: 0)) if as_dict else defaultdict(list)
        if as_dict:
            def putz(z, x):
                locs[z][x] += 1
        else:
            def putz(z, x):
                locs[z].append(x)

        G = json.load(open(path,'r'))
        if len(G['Isolate_names']) < 10:
            continue

        for db in G['Nodes'].values():
            gs = [db['Genomes']['Consensus']] if consensus else db['Genomes']['Alignment'].values()
            for g in gs:
                Z_fwd, Z_rev = 0, 0

                # Build first kmer 
                nshift = 0
                for i in range(k):
                    c = nuctab[ord(g[i])]
                    if c > 4:
                        nshift = 0
                        continue
                    Z_fwd, Z_rev = shift(c, Z_fwd, Z_rev)
                    nshift += 1

                if nshift == k:
                    putz(Z_fwd, 0)
                    putz(Z_rev, 0)

                # Build successive kmers by iteratively shifting
                nshift = 0
                L = len(g)
                for i in range(k, L):
                    c = nuctab[ord(g[i])]
                    if c > 4:
                        nshift = 0
                        continue

                    Z_fwd, Z_rev = shift(c, Z_fwd, Z_rev)
                    nshift += 1
                    if nshift >= k:
                        l = 1 - abs(1 - 2*(i-k+1)/L) if normalize else L/2 - abs(L/2 - (i - k + 1))
                        putz(Z_fwd, l)
                        putz(Z_rev, l)

        kmer_locs[path] = default_to_normal(locs)
        n += 1
        if n >= 30:
           break

    return kmer_locs

def to_kmer(Z, k):
    # Elementwise function
    def xform(z):
        inv  = np.array(['A', 'C', 'G', 'T'], dtype=np.str)
        kmer = np.empty(k, dtype=np.str)
        mask = 3
        for i in range(k):
            kmer[k-i-1] = inv[mask & z]
            z >>= 2
        return kmer

    if isinstance(Z, Iterable):
        return [xform(z) for z in Z]
    else:
        return xform(Z)

def from_kmer(kmer):
    assert len(kmer) < 32

    Z_fwd, Z_rev = 0, 0
    shift  = shiftfn(len(kmer))
    nshift = 0
    for nuc in kmer:
        c = nuctab[ord(nuc)]
        if c > 4:
            nshift = 0
            continue
        Z_fwd, Z_rev = shift(c, Z_fwd, Z_rev)
        nshift += 1

    if nshift == len(kmer):
        return Z_fwd, Z_rev

    return None, None


# ------------------------------------------------------------------------
# Kmer at edges functions

def kmers_at_boundary(k=8, delta=500):
    assert k < delta

    edge_kmers = defaultdict(lambda: 0)
    midl_kmers = defaultdict(lambda: 0)
    shift = shiftfn(k)

    def sketch(kmers, seq, lb, ub):
        nshift = 0
        Z_fwd, Z_rev = 0, 0

        for i in range(lb, lb+k):
            c = nuctab[ord(seq[i])]
            if c > 4:
                nshift = 0
                continue
            Z_fwd, Z_rev = shift(c, Z_fwd, Z_rev)
            nshift += 1

        if nshift == k:
            kmers[frozenset({Z_fwd, Z_rev})] += 1

        # Build successive kmers by iteratively shifting
        for i in range(lb+k, ub):
            c = nuctab[ord(g[i])]
            if c > 4:
                nshift = 0
                continue

            Z_fwd, Z_rev = shift(c, Z_fwd, Z_rev)
            nshift += 1
            if nshift >= k:
                kmers[frozenset({Z_fwd, Z_rev})] += 1

        return kmers

    for path in glob(f"{visdir}/*.json"):
        G = json.load(open(path,'r'))
        if len(G['Isolate_names']) < 10:
            continue

        print(f"----> Analyzing {path}")

        for db in G['Nodes'].values():
            g = db['Genomes']['Consensus']
            L = len(g)
            if L < 5*delta:
                continue

            edge_kmers = sketch(edge_kmers, g, 0, delta)
            edge_kmers = sketch(edge_kmers, g, L-delta, L)
            midl_kmers = sketch(midl_kmers, g, L//2-delta, L//2)
            midl_kmers = sketch(midl_kmers, g, L//2, L//2 + delta)

    # Build first kmer 
    return edge_kmers, midl_kmers

# ------------------------------------------------------------------------
# Start codon position

def start_codon_pos(normalize=True):
    codon_pos = [[], [], []]

    for path in glob(f"{visdir}/*.json"):
        print(f"----> Analyzing {path}")

        G = json.load(open(path,'r'))
        if len(G['Isolate_names']) < 10:
            continue

        for db in G['Nodes'].values():
            g = db['Genomes']['Consensus']

            L = len(g)
            for i in range(3):
                for x in range(i, L, 3):
                    if g[x:x+3] == "ATG":
                        l = 1 - abs(1 - 2*(x-3+1)/L) if normalize else L/2 - abs(L/2 - (x-3+1))
                        codon_pos[i].append(l)

    return codon_pos

# ------------------------------------------------------------------------
# Utility plotting functions

def cdfplot(x, **kwargs):
    plt.plot(sorted(x), np.linspace(0, 1, len(x)), **kwargs)

# -----------------------------------------------------------------------------
# Main point of entry

K = 10
D = 400
if __name__ == "__main__":
### Get clonality and length of blocks
    # divs = clonality_and_length()

### Bin all kmers based on fractional distance from block ends.
    print("--> Getting kmer positions")
    kmer_locs = bin_kmer_location(k=K, normalize=False)
    print("--> Merging dictionaries")
    all_locs = merge_dicts(kmer_locs)

    print("--> Computing medians")
    data = [ (np.median(x), len(x), Z) for Z, x in all_locs.items() ]
    M, N, Z = zip(*data)
    M = np.array(M)
    N = np.array(N)
    Z = np.array(Z)

    plt.scatter(M, N, alpha=.05)
    plt.yscale('log')
    plt.xlabel("Median kmer position from end")
    plt.ylabel("Number of occurences")
    plt.title(f"k = {K}")
    cm = mpl.cm.get_cmap("Spectral")
    for x in all_locs.values():
        if len(x) < 20:
            continue
        plt.plot(sorted(x), np.linspace(0, 1, len(x)), color=cm(np.median(x)), alpha=.2)

### Bin all kmers only delta bp away from ends
    # kmers_e, kmers_m = kmers_at_boundary(k=K, delta=D)

    # cdfplot(kmers_e.values(), label=f"{D} bp from block edge")
    # cdfplot(kmers_m.values(), label=f"{D} bp from block middle")
    # plt.xlabel("Number of occurrences")
    # plt.xscale('log')
    # plt.ylabel("CDF")
    # plt.title(f"k={K}")
    # plt.legend()

### Get start codon positions
    # start_pos = start_codon_pos(normalize=False)
