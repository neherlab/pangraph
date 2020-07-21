import os, sys
import numpy as np

from glob  import glob

from Bio           import SeqIO, Phylo
from Bio.Seq       import Seq
from Bio.SeqRecord import SeqRecord

from .         import suffix
from .block    import Block
from .sequence import Node, Path
from .utils    import Strand, as_string, parse_paf, panic, tryprint, asrecord, new_strand, breakpoint

# ------------------------------------------------------------------------
# Global variables

outdir      = "data/graph"
MAXSELFMAPS = 100

# ------------------------------------------------------------------------
# Graph class

class Graph(object):
    """docstring for Graph"""

    def __init__(self):
        self.name = ""   # The name of graph. Will be used as file basename in exports
        self.blks = {}   # All blocks/alignments
        self.seqs = {}   # All sequences (as list of blocks)
        self.sfxt = None # Suffix tree of block records
            def from_dict():
        self.dmtx = None # Graph distance matrix

    # --- Class methods ---

    @classmethod
    def from_seq(cls, name, seq):
        newg = cls()
        blk  = Block.from_seq(name, seq)
        newg.name = name
        newg.blks = {blk.id : blk}
        newg.seqs = {name : Path((blk.id, Strand.Plus, 0), 0)}

        return newg

    @classmethod
    def from_dict(cls, d):
        G = Graph()
        G.name = d['name']
        G.blks = [Block.from_dict(b) for b in d['blocks']]
        G.blks = {b.id : b for b in G.blks}
        G.seqs = {key:Path.from_dict(val) for key,val in d['seqs'].items()}
        G.sfxt = None
        G.dmtx = None
        if d['suffix'] is not None:
            G.compile_suffix()
            G.dmtx = d['distmtx']

        return G

    @classmethod
    def connected_components(cls, G):
        # internal functions
        def overlaps(s1, s2):
            return len(s1.intersection(s2)) > 0
        def component(graph, name):
            cc = Graph()
            cc.blks = {id:G.blks.pop(id) for id in graph}
            cc.seqs = {nm:G.seqs.pop(nm) for nm in name}
            cc.sfxt = None
            cc.dmtx = None
            return cc

        # main body
        graphs, names = [], []
        for name, seq in G.seqs.items():
            blks = set(s[0] for s in seq)
            gi   = [ i for i, g in enumerate(graphs) if overlaps(blks, g)]
            if len(gi) == 0:
                graphs.append(blks)
                names.append(set([name]))
                continue

            graphs[gi[0]] = graphs[gi[0]].union(blks, *(graphs.pop(i) for i in gi[:0:-1]))
            names[gi[0]]  = names[gi[0]].union(set([name]), *(names.pop(i) for i in gi[:0:-1]))

        return [component(graph, name) for graph, name in zip(graphs, names)]

    @classmethod
    def fuse(cls, g1, g2):
        ng = Graph()
        ng.blks = {}
        ng.blks.update(g1.blks)
        ng.blks.update(g2.blks)

        ng.seqs = {s:list(b) for s,b in list(g1.seqs.items())+list(g2.seqs.items())}

        return ng

    # ---------------
    # methods

    def union(self, qpath, rpath, out, cutoff=0, mu=10, beta=2, extensive=False):
        import warnings
        from skbio.alignment import global_pairwise_align
        from skbio import DNA

        # ----------------------------------
        # internal functions

        def global_aln(s1, s2):
            M = {}
            alphas = ['A', 'C', 'G', 'T', 'N']
            for a in alphas:
                M[a] = {}
                for b in alphas:
                    if a == b:
                        M[a][b] = 2
                    else:
                        M[a][b] = -3
            return global_pairwise_align(s1, s2, 5, 2, M)

        def energy(hit):
            l    = hit["aligned_bases"]
            if l <= cutoff:
                return l

            num  = lambda k: len(self.blks[hit[k]["name"]].muts)
            cuts = lambda k: (hit[k]['start'] > cutoff) + ((hit[k]['len']-hit[k]['end']) > cutoff)

            if extensive:
                delP = num('qry')*cuts('qry') + num('ref')*cuts('ref')
            else:
                delP = cuts('qry') + cuts('ref')
            dmut = hit["aligned_length"] * hit["divergence"]

            return -l + mu*delP + beta*dmut

        def accepted(hit):
            return energy(hit) < 0

        if cutoff <= 0:
            def proc(hit):
                return hit
        else:
            def proc(hit):
                # -----------------------
                # load in sequences
                with open(f"{qpath}.fa", 'r') as fd:
                    qfa = {s.id:str(s.seq) for s in SeqIO.parse(fd, 'fasta')}

                if qpath == rpath:
                    rfa = qfa
                else:
                    with open(f"{rpath}.fa", 'r') as fd:
                        rfa = {s.id:str(s.seq) for s in SeqIO.parse(fd, 'fasta')}

                # -----------------------
                # internal functions

                def tocigar(aln):
                    cigar = ""
                    s1, s2 = np.fromstring(aln[0], dtype=np.int8), np.fromstring(aln[1], dtype=np.int8)
                    M, I, D = 0, 0, 0
                    for (c1, c2) in zip(s1, s2):
                        if c1 == ord("-") and c2 == ord("-"):
                            print("PANIC")
                            import ipdb; ipdb.set_trace()
                        elif c1 == ord("-"):
                            if I > 0:
                                cigar += f"{I}I"
                                I = 0
                            elif M > 0:
                                cigar += f"{M}M"
                                M = 0
                            D += 1
                        elif c2 == ord("-"):
                            if D > 0:
                                cigar += f"{D}D"
                                D = 0
                            elif M > 0:
                                cigar += f"{M}M"
                                M = 0
                            I += 1
                        else:
                            if D > 0:
                                cigar += f"{D}D"
                                D = 0
                            elif I > 0:
                                cigar += f"{I}I"
                                I = 0
                            M += 1
                    if I > 0:
                        cigar += f"{I}I"
                        I = 0
                    elif M > 0:
                        cigar += f"{M}M"
                        M = 0
                    elif D > 0:
                        cigar += f"{D}D"
                        M = 0

                    return cigar

                def revcmpl_if(s, cond):
                    if cond:
                        return str(Seq.reverse_complement(Seq(s)))
                    else:
                        return s

                def getseqs():
                    return qfa[hit['qry']['name']], rfa[hit['ref']['name']]

                # -----------------------
                # body

                dS_q = hit['qry']['start']
                dE_q = hit['qry']['len'] - hit['qry']['end']
                dS_r = hit['ref']['start']
                dE_r = hit['ref']['len'] - hit['ref']['end']

                warnings.simplefilter("ignore")

                # if hit['qry']['start'] > 0 and hit['qry']['start'] < 50 or \
                #    hit['ref']['start'] > 0 and hit['ref']['start'] < 50:
                #     breakpoint("START CHANGE")

                # Left side of match
                if 0 < dS_q <= cutoff and (dS_r > cutoff or dS_r == 0):
                    hit['cigar'] = f"{dS_q}I" + hit['cigar']
                    hit['qry']['start'] = 0
                elif 0 < dS_r <= cutoff and (dS_q > cutoff or dS_q == 0):
                    hit['cigar'] = f"{dS_r}D" + hit['cigar']
                    hit['ref']['start'] = 0
                elif 0 < dS_q <= cutoff and 0 < dS_r <= cutoff:
                    qseq, rseq = getseqs()
                    aln = global_aln(DNA(revcmpl_if(qseq, hit['orientation']==Strand.Minus)[0:dS_q]), DNA(rseq[0:dS_r]))
                    aln = tuple([str(v) for v in aln[0].to_dict().values()])
                    hit['cigar'] = tocigar(aln) + hit['cigar']
                    hit['qry']['start'] = 0
                    hit['ref']['start'] = 0
                    hit['aligned_bases'] += len(aln[0])

                # Right side of match
                if 0 < dE_q <= cutoff and (dE_r > cutoff or dE_r == 0):
                    hit['cigar'] += f"{dE_q}I"
                    hit['qry']['end'] = hit['qry']['len']
                elif 0 < dE_r <= cutoff and (dE_q > cutoff or dE_q == 0):
                    hit['cigar'] += f"{dE_r}D"
                    hit['ref']['end'] = hit['ref']['len']
                elif 0 < dE_q <= cutoff and 0 < dE_r <= cutoff:
                    qseq, rseq = getseqs()
                    aln = global_aln(DNA(revcmpl_if(qseq, hit['orientation']==Strand.Minus)[-dE_q:]), DNA(rseq[-dE_r:]))
                    aln = tuple([str(v) for v in aln[0].to_dict().values()])

                    hit['cigar'] = hit['cigar'] + tocigar(aln)
                    hit['qry']['end'] = hit['qry']['len']
                    hit['ref']['end'] = hit['ref']['len']
                    hit['aligned_bases'] += len(aln[0])

                return hit

        # ----------------------------------
        # body

        os.system(f"minimap2 -t 2 -x asm20 -m 10 -n 2 -s 30 -D -c {rpath}.fa {qpath}.fa 1>{out}.paf 2>log")

        with open(f"{out}.paf") as fd:
            paf = parse_paf(fd)
        paf.sort(key = lambda x: energy(x))

        merged_blks = set()
        if len(paf) == 0:
            return self, False

        merged = False
        for hit in paf:
            if hit['qry']['name'] in merged_blks \
            or hit['ref']['name'] in merged_blks \
            or(hit['ref']['name'] <= hit['qry']['name'] and qpath == rpath) \
            or not accepted(hit):
                continue

            tryprint(f"------> merge {hit['ref']['name']} with {hit['qry']['name']}", False)
            self.merge(proc(hit))
            merged = True
            merged_blks.add(hit['ref']['name'])
            merged_blks.add(hit['qry']['name'])

        self.purge_empty()
        return self, merged

    def prune(self):
        blks_remain = set()
        for iso in self.seqs:
            blks_remain.update([s[0] for s in self.seqs[iso]])
        self.blks = {b:self.blks[b] for b in blks_remain}

        return

    def purge_empty(self):
        for iso in self.seqs:
            goodblks = []
            popblks  = set()
            for i, (b, _, n) in enumerate(self.seqs[iso]):
                if b in popblks:
                    continue

                if self.blks[b].is_empty(iso, n):
                    if (iso, n) not in self.blks[b].muts:
                        import ipdb; ipdb.set_trace()
                    self.blks[b].muts.pop((iso, n))
                else:
                    goodblks.append(i)

                if not self.blks[b].has(iso):
                    popblks.add(b)

            self.seqs[iso] = [self.seqs[iso][i] for i in goodblks]

        return

    def merge(self, hit):
        refblk_orig = self.blks[hit['ref']['name']]
        qryblk_orig = self.blks[hit['qry']['name']]

        # As we slice here, we DONT need to remember the starting position.
        # This is why in from_aln(aln) we set the start index to 0
        refblk = refblk_orig[hit['ref']['start']:hit['ref']['end']]
        qryblk = qryblk_orig[hit['qry']['start']:hit['qry']['end']]

        if hit["orientation"] == Strand.Minus:
            qryblk = qryblk.rev_cmpl()

        aln = {"ref_seq"     : as_string(refblk.seq), # "".join(refblk.seq),
               "qry_seq"     : as_string(qryblk.seq), # "".join(qryblk.seq),
               "cigar"       : hit["cigar"],
               "ref_cluster" : refblk.muts,
               "qry_cluster" : qryblk.muts,
               "ref_start"   : hit["ref"]["start"],
               "ref_name"    : hit["ref"]["name"],
               "qry_start"   : hit["qry"]["start"],
               "qry_name"    : hit["qry"]["name"],
               "orientation" : hit["orientation"]}

        merged_blks, qryblks, refblks, isomap = Block.from_aln(aln)
        for merged_blk in merged_blks:
            self.blks[merged_blk.id] = merged_blk

        def update(blk, addblks, hit, strand):
            new_blks = []

            # The convention for the tuples are (block id, strand orientation, if merged)
            if hit['start'] > 0:
                if hit['start'] < 50:
                    breakpoint("TOO SMALL")

                left = blk[0:hit['start']]
                self.blks[left.id] = left
                new_blks.append((left.id, Strand.Plus, False))

            for b in addblks:
                new_blks.append((b.id, strand, True))

            if hit['end'] < len(blk):
                right = blk[hit['end']:]
                self.blks[right.id] = right
                new_blks.append((right.id, Strand.Plus, False))

            def replace(tag, blk, newblks):
                iso = tag[0]

                new_blk_seq = []
                oldseq = self.seqs[iso].copy()

                for b in self.seqs[iso]:
                    if b[0] == blk.id and b[2] == tag[1]:
                        orig_strand    = b[1]
                        tmp_new_blocks = [(ID, new_strand(orig_strand, ns), isomap[ID][blk.id][tag][1]) if merged else (ID, new_strand(orig_strand, ns), b[2]) for ID, ns, merged in newblks]
                        if orig_strand == Strand.Minus:
                            tmp_new_blocks = tmp_new_blocks[::-1]

                        new_blk_seq.extend(tmp_new_blocks)
                    else:
                        new_blk_seq.append(b)

                self.seqs[iso] = new_blk_seq
                for b in self.seqs[iso]:
                    if (iso, b[2]) not in self.blks[b[0]].muts:
                        import ipdb; ipdb.set_trace()


            for tag in blk.muts.keys():
                replace(tag, blk, new_blks)

        update(refblk_orig, refblks, hit['ref'], Strand.Plus)
        update(qryblk_orig, qryblks, hit['qry'], hit['orientation'])
        self.prune()

        return

    def extract(self, name, strip_gaps=True, verbose=False):
        seq = ""
        for (b, strand, num) in self.seqs[name]:
            tmp_seq = self.blks[b].extract(name, num, strip_gaps=False, verbose=verbose)
            if strand == Strand.Plus:
                seq += tmp_seq
            else:
                seq += str(Seq.reverse_complement(Seq(tmp_seq)))

        # start_pos = self.spos[name]

        if start_pos:
            seq = seq[start_pos:] + seq[:start_pos]

        if strip_gaps:
            seq = seq.replace('-', '')

        return seq

    def compress_ratio(self, extensive=False, name=None):
        unc = 0
        if name is None:
            for n in self.seqs:
                seq  = self.extract(n)
                unc += len(seq)
            cmp = np.sum([len(x) for x in self.blks.values()])
        else:
            cmp = np.sum([len(x) for x in self.blks.values() if name in x.muts])

        return unc/cmp/len(self.seqs) if not extensive else unc/cmp

    def contains(self, other):
        return set(other.seqs.keys()).issubset(set(self.seqs.keys()))

    def compile_suffix(self, force=False):
        if self.sfxt is None or force:
            self.sfxt = suffix.Tree({k: [c[0:2] for c in v] for k, v in self.seqs.items()})

    def compute_pdist(self, force=False):
        if self.dmtx is None or force:
            nms, N = sorted(list(self.seqs.keys())), len(self.seqs)
            self.dmtx = np.zeros((N*(N-1))//2)

            n = 0
            for i, nm1 in enumerate(nms):
                for nm2 in nms[:i]:
                    self.dmtx[n] = len(self.sfxt.matches(nm1, nm2))
                    n += 1

    def to_json(self, wtr, minlen=500):
        import json

        J = {}
        cleaned_seqs = {s:[b for b in self.seqs[s] if len(self.blks[b[0]])>minlen]
                             for s in self.seqs}
        relevant_blocks = set()
        for s in cleaned_seqs.values():
            relevant_blocks.update([b[0] for b in s])

        J['Isolate_names'] = list(cleaned_seqs.keys())
        J['Plasmids']      = [[x for x in cleaned_seqs[s]] for s in J['Isolate_names']]

        # Build node substructure
        nodes = {}
        for b in relevant_blocks:
            aln = { J["Isolate_names"].index(iso) :
                    self.blks[b].extract(iso, num, strip_gaps=False) for iso, num in self.blks[b].muts }
            nodes[b] = {"ID"        : b,
                        "Genomes"   : {"Consensus" : ''.join(self.blks[b].seq),
                                       "Alignment" : aln },
                        "Out_Edges" : [],
                        "In_Edges"  : []}

        # Gather edges (node/node junctions) for each node
        edges = {}
        for pname, p in zip(range(len(J["Isolate_names"])), J['Plasmids']):
            for i in range(len(p)-1):
                e = (p[i], p[i+1])
                if e in edges:
                    edges[e]["Isolates"].append(pname)
                else:
                    edges[e] = {"Source" : e[0], "Target" : e[1], "Isolates" : [pname]}
            e = (p[-1], p[0])
            if e in edges:
                edges[e]["Isolates"].append(pname)
            else:
                edges[e] = {"Source" : e[0], "Target" : e[1], "Isolates" : [pname]}

        for e in edges:
            nodes[e[0][0]]["Out_Edges"].append(edges[e])
            nodes[e[1][0]]["In_Edges"].append(edges[e])

        J["Nodes"] = nodes
        json.dump(J, wtr)

        return

    def to_dict(self):
        return {'name'   : self.name,
                'seqs'   : self.seqs,
                'blocks' : [b.to_dict() for b in self.blks.values()],
                'suffix' : None if self.sfxt is None else "compiled",
                'distmtx': self.dmtx}

    def write_fasta(self, wtr):
        SeqIO.write(sorted([ SeqRecord(seq=Seq(as_string(c.seq)), id=c.id, description='')
            for c in self.blks.values() ], key=lambda x: len(x), reverse=True), wtr, format='fasta')

        return

