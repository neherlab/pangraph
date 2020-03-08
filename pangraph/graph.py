import os, sys
import numpy as np
import pyfaidx as fai

from glob  import glob

from Bio           import SeqIO, Phylo
from Bio.Seq       import Seq
from Bio.SeqRecord import SeqRecord

from . import suffix
from .block import Block
from .utils import Strand, asstring, parsepaf, panic, tryprint, asrecord, newstrand

# ------------------------------------------------------------------------
# Global variables
outdir      = "data/graph"
maxselfmaps = 100

# ------------------------------------------------------------------------
# Graph class

class Graph(object):
    """docstring for Graph"""

    visdir = f"{outdir}/vis"
    alndir = f"{outdir}/aln"
    blddir = f"{outdir}/bld"

    def __init__(self):
        super(Graph, self).__init__()
        self.name = ""   # The name of graph. Will be used as file basename in exports
        self.blks = {}   # All blocks/alignments
        self.seqs = {}   # All sequences (as list of blocks)
        self.spos = {}   # All start positions (index of block that starts fasta)
        self.sfxt = None # Suffix tree of block records
        self.dmtx = None # Graph distance matrix

    # --- Class methods ---

    @classmethod
    def fromseq(cls, name, seq):
        newg = cls()
        blk  = Block.fromseq(name, seq)
        newg.name = name
        newg.blks = {blk.id : blk}
        newg.seqs = {name : [(blk.id, Strand.Plus, 0)]}
        newg.spos = {name : 0}

        return newg

    @classmethod
    def fromdict(cls, d):
        G = Graph()
        G.name = d['name']
        G.blks = [Block.fromdict(b) for b in d['blocks']]
        G.blks = {b.id : b for b in G.blks}
        G.seqs = d['seqs']
        G.spos = d['starts']
        G.sfxt = None
        G.dmtx = None
        if d['suffix'] is not None:
            G.compilesuffix()
            G.dmtx = d['distmtx']

        return G

    @classmethod
    def fuse(cls, g1, g2):
        ng = Graph()
        ng.blks = {}
        ng.blks.update(g1.blks)
        ng.blks.update(g2.blks)

        ng.seqs = {s:list(b) for s,b in list(g1.seqs.items())+list(g2.seqs.items())}
        ng.spos = {s:b for s,b in list(g1.spos.items())+list(g2.spos.items())}

        return ng

    @classmethod
    def fromnwk(cls, path, seqs, save=True, verbose=False, mu=100, beta=2):
        # Debugging function that will check reconstructed sequence against known real one.
        def check(seqs, T, G, verbose=False):
            nerror = 0
            uncompressed_length = 0
            for n in T.get_terminals():
                if n.name not in G.seqs:
                    continue

                seq = seqs[n.name]
                orig = str(seq[:].seq).upper()
                tryprint(f"--> Checking {n.name}", verbose=verbose)
                rec = G.extract(n.name)
                uncompressed_length += len(orig)
                if orig != rec:
                    nerror += 1

                    with open("test.fa", "w+") as out:
                        out.write(f">original\n{orig}\n")
                        out.write(f">reconstructed\n{rec}")

                    for i in range(len(orig)//100):
                        if (orig[i*100:(i+1)*100] != rec[i*100:(i+1)*100]):
                            print("-----------------")
                            print("O:", i, orig[i*100:(i+1)*100])
                            print("G:", i, rec[i*100:(i+1)*100])

                            diffs = [i for i in range(len(rec)) if rec[i] != orig[i]]
                            pos   = [0]
                            blks  = G.seqs[n.name]
                            for b, strand, num in blks:
                                pos.append(pos[-1] + len(G.blks[b].extract(n.name, num)))
                            pos = pos[1:]

                            testseqs = []
                            for b in G.seqs[n.name]:
                                if b[1] == Strand.Plus:
                                    testseqs.append("".join(G.blks[b[0]].extract(n.name, b[2])))
                                else:
                                    testseqs.append("".join(Seq.reverse_complement(G.blks[b[0]].extract(n.name, b[2]))))

                            import ipdb; ipdb.set_trace()
                else:
                    tryprint(f"+++ Verified {n.name}", verbose=verbose)

            if nerror == 0:
                tryprint("all sequences correctly reconstructed", verbose=verbose)
                tlength = np.sum([len(x) for x in G.blks.values()])
                tryprint(f"--- total graph length: {tlength}", verbose=verbose)
                tryprint(f"--- total input sequence: {uncompressed_length}", verbose=verbose)
                tryprint(f"--- compression: {uncompressed_length/tlength:1.2f}", verbose=verbose)
            else:
                raise ValueError("bad sequence reconstruction")

        T = Phylo.read(path, "newick")
        if T.count_terminals() == 1:
            return Graph(), False

        for i, n in enumerate(T.get_terminals()):
            seq      = asrecord(seqs[n.name])
            n.graph  = Graph.fromseq(seq.id, str(seq.seq).upper())
            n.name   = seq.id
            n.fapath = f"{Graph.blddir}/{n.name}"
            tryprint(f"------> Outputting {n.fapath}", verbose=verbose)
            n.graph.tofasta(n.fapath)

        nnodes = 0
        for n in T.get_nonterminals(order='postorder'):
            nnodes += 1

            # Simple graph "fuse". Straight concatenation
            n.name   = f"NODE_{nnodes:07d}"
            n.graph  = Graph.fuse(n.clades[0].graph, n.clades[1].graph)
            # check(seqs, T, n.graph)
            n.fapath = os.path.join(*[Graph.blddir, f"{n.name}.fasta"])

            tryprint(f"-- node {n.name} with {len(n.clades)} children", verbose)

            # Initial map from graph to itself
            n.graph, _ = n.graph._mapandmerge(n.clades[0].fapath, n.clades[1].fapath, f"{Graph.blddir}/{n.name}",
                            cutoff=50, mu=mu, beta=beta)

            i, contin = 0, True
            while contin:
                tryprint(f"----> merge round {i}", verbose)
                # check(seqs, T, n.graph)
                itrseq = f"{Graph.blddir}/{n.name}_iter_{i}"
                n.graph.tofasta(itrseq)
                n.graph, contin = n.graph._mapandmerge(itrseq, itrseq, f"{Graph.blddir}/{n.name}_iter_{i}",
                            cutoff=50, mu=mu, beta=beta)
                i += 1

                contin &= i < maxselfmaps

            tryprint(f"-- Blocks: {len(n.graph.blks)}, length: {np.sum([len(b) for b in n.graph.blks.values()])}\n", verbose)
            n.graph.tofasta(f"{Graph.blddir}/{n.name}")
            # Continuous error logging
            # print(f"Node {n.name}", file=sys.stderr)
            # print(f"--> Compression ratio parent: {n.graph.compressratio()}", file=sys.stderr)
            # print(f"--> Compression ratio child1: {n.clades[0].graph.compressratio()}", file=sys.stderr)
            # print(f"--> Compression ratio child2: {n.clades[1].graph.compressratio()}", file=sys.stderr)
            # Output with content
            # print(f"{n.graph.compressratio()}\t{n.clades[0].graph.compressratio()}\t{n.clades[1].graph.compressratio()}\t{n.clades[0].branch_length + n.clades[1].branch_length}", file=sys.stdout, flush=True)

        G      = T.root.graph
        G.name = ".".join(os.path.basename(path).split(".")[:-1])
        # check(seqs, T, G)

        if save:
            G.tojson()
            G.tofasta()

        return G, True

    @classmethod
    def cleanbld(cls):
        for p in glob(f"{cls.blddir}/*"):
            os.remove(p)

    # --- Internal methods ---

    def _mapandmerge(self, qpath, rpath, out, cutoff=None, mu=100, beta=2):
        # from Bio import pairwise2
        import warnings
        from skbio.alignment import global_pairwise_align
        from skbio import DNA

        def globalaln(s1, s2):
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
            delP = len(self.blks[hit["qry"]["name"]].muts) + \
                   len(self.blks[hit["ref"]["name"]].muts)

            dmut = hit["aligned_length"] * hit["divergence"]

            # print(f"{l} :: {mu*delP + beta*dmut}")

            return -l + mu*delP + beta*dmut

        def accepted(hit):
            return energy(hit) < 0

        qpath = qpath.replace(".fasta", "")
        rpath = rpath.replace(".fasta", "")

        # print(f"--> {qpath} & {rpath} -> {out}.paf")
        os.system(f"minimap2 -t 2 -x asm5 -D -c {rpath}.fasta {qpath}.fasta 1>{out}.paf 2>log")

        paf = parsepaf(f"{out}.paf")
        paf.sort(key = lambda x: energy(x))

        merged_blks = set()
        if len(paf) == 0:
            return self, False

        if cutoff is None:
            def proc(hit):
                return hit
        else:
            def proc(hit):
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
                    if qpath == rpath:
                        fa   = fai.Fasta(f"{qpath}.fasta")
                        qseq = fa[hit['qry']['name']][:].seq
                        rseq = fa[hit['ref']['name']][:].seq
                    else:
                        qseq = fai.Fasta(f"{qpath}.fasta")[hit['qry']['name']][:].seq
                        rseq = fai.Fasta(f"{rpath}.fasta")[hit['ref']['name']][:].seq

                    return qseq, rseq

                dS_q = hit['qry']['start']
                dE_q = hit['qry']['len'] - hit['qry']['end']
                dS_r = hit['ref']['start']
                dE_r = hit['ref']['len'] - hit['ref']['end']

                warnings.simplefilter("ignore")

                # Left side of match
                if 0 < dS_q <= cutoff and dS_r > cutoff:
                    hit['cigar'] = f"{dS_q}I" + hit['cigar']
                    hit['qry']['start'] = 0
                elif 0 < dS_r <= cutoff and dS_q > cutoff:
                    hit['cigar'] = f"{dS_r}D" + hit['cigar']
                    hit['ref']['start'] = 0
                elif 0 < dS_q <= cutoff and 0 < dS_r <= cutoff:
                    qseq, rseq = getseqs()
                    aln = globalaln(DNA(revcmpl_if(qseq, hit['orientation']==Strand.Minus)[0:dS_q]), DNA(rseq[0:dS_r]))
                    aln = tuple([str(v) for v in aln[0].to_dict().values()])
                    hit['cigar'] = tocigar(aln) + hit['cigar']
                    hit['qry']['start'] = 0
                    hit['ref']['start'] = 0
                    hit['aligned_bases'] += len(aln[0])

                # Right side of match
                if 0 < dE_q <= cutoff and dE_r > cutoff:
                    hit['cigar'] += f"{dE_q}I"
                    hit['qry']['end'] = hit['qry']['len']
                elif 0 < dE_r <= cutoff and dE_q > cutoff:
                    hit['cigar'] += f"{dE_r}D"
                    hit['ref']['end'] = hit['ref']['len']
                elif 0 < dE_q <= cutoff and 0 < dE_r <= cutoff:
                    qseq, rseq = getseqs()
                    aln = globalaln(DNA(revcmpl_if(qseq, hit['orientation']==Strand.Minus)[-dE_q:]), DNA(rseq[-dE_r:]))
                    aln = tuple([str(v) for v in aln[0].to_dict().values()])

                    hit['cigar'] = hit['cigar'] + tocigar(aln)
                    hit['qry']['end'] = hit['qry']['len']
                    hit['ref']['end'] = hit['ref']['len']
                    hit['aligned_bases'] += len(aln[0])

                return hit

        merged = False
        for hit in paf:
            if hit['qry']['name'] in merged_blks \
            or hit['ref']['name'] in merged_blks \
            or(hit['ref']['name'] <= hit['qry']['name'] and qpath == rpath) \
            or hit['align_score'] < 0.5:
                continue

            if not accepted(hit):
                continue

            tryprint(f"------> merge {hit['ref']['name']} with {hit['qry']['name']}", False)
            self.merge(proc(hit))
            merged = True
            merged_blks.add(hit['ref']['name'])
            merged_blks.add(hit['qry']['name'])

        self.purgeempty()
        return self, merged

    # --- Instance methods ---

    def prune(self):
        blks_remain = set()
        for iso in self.seqs:
            blks_remain.update([s[0] for s in self.seqs[iso]])
        self.blks = {b:self.blks[b] for b in blks_remain}

        return

    def purgeempty(self):
        for iso in self.seqs:
            goodblks = []
            popblks  = set()
            for i, (b, _, n) in enumerate(self.seqs[iso]):
                if b in popblks:
                    continue

                # NOTE: Concatenating the sequence just to see if it's empty
                #       takes 1/2 the time of running! 

                # bseq = self.blks[b].extract(iso, n)
                # if bseq:
                #     goodblks.append(i)
                # else:
                #     self.blks[b].muts.pop((iso, n))

                if self.blks[b].isempty(iso, n):
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
            qryblk = qryblk.revcmpl()

        aln = {"ref_seq"     : asstring(refblk.seq), # "".join(refblk.seq),
               "qry_seq"     : asstring(qryblk.seq), # "".join(qryblk.seq),
               "cigar"       : hit["cigar"],
               "ref_cluster" : refblk.muts,
               "qry_cluster" : qryblk.muts,
               "ref_start"   : hit["ref"]["start"],
               "ref_name"    : hit["ref"]["name"],
               "qry_start"   : hit["qry"]["start"],
               "qry_name"    : hit["qry"]["name"],
               "orientation" : hit["orientation"]}

        # shares_isos = not set(self.blks[aln['ref_name']].muts.keys()).isdisjoint(set(self.blks[aln['qry_name']].muts.keys()))
        merged_blks, qryblks, refblks, isomap = Block.fromaln(aln) #, debug=shares_isos)
        for merged_blk in merged_blks:
            self.blks[merged_blk.id] = merged_blk

        def update(blk, addblks, hit, strand):
            new_blks = []
            # The convention for the tuples are (block id, strand orientation, if merged)

            if hit['start'] > 0:
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
                        tmp_new_blocks = [(ID, newstrand(orig_strand, ns), isomap[ID][blk.id][tag][1]) if merged else (ID, newstrand(orig_strand, ns), b[2]) for ID, ns, merged in newblks]
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

        start_pos = self.spos[name]

        if start_pos:
            seq = seq[start_pos:] + seq[:start_pos]

        if strip_gaps:
            seq = seq.replace('-', '')

        return seq

    def compressratio(self, name=None):
        unclen = 0
        if name is None:
            for n in self.seqs:
                seq = self.extract(n)
                unclen += len(seq)

        if name is None:
            cmplen = np.sum([len(x) for x in self.blks.values()])
        else:
            cmplen = np.sum([len(x) for x in self.blks.values() if name in x.muts])

        return unclen/cmplen

    def compilesuffix(self, force=False):
        if self.sfxt is None or force:
            self.sfxt = suffix.Tree({k: [c[0:2] for c in v] for k, v in self.seqs.items()})

    def computepairdists(self, force=False):
        if self.dmtx is None or force:
            nms, N = sorted(list(self.seqs.keys())), len(self.seqs)
            self.dmtx = np.zeros((N*(N-1))//2)

            n = 0
            for i, nm1 in enumerate(nms):
                for nm2 in nms[:i]:
                    self.dmtx[n] = len(self.sfxt.matches(nm1, nm2))
                    n += 1

    def tojson(self, minlen=500):
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

        # print(f"----> Saving to {Graph.visdir}/{self.name}.json")
        with open(f"{Graph.visdir}/{self.name}.json", 'w+') as fh:
            json.dump(J, fh)

        return

    def tofasta(self, path=None):
        if path is None:
            path = f"{Graph.alndir}/{self.name}.fasta"
        elif not path.endswith(".fasta"):
            path = path + ".fasta"
        SeqIO.write(sorted([ SeqRecord(seq=Seq(asstring(c.seq)), id=c.id, description='')
            for c in self.blks.values() ], key=lambda x: len(x), reverse=True), path, format='fasta')

        return

    def todict(self):
        return {'name'   : self.name,
                'seqs'   : self.seqs,
                'starts' : self.spos,
                'blocks' : [b.todict() for b in self.blks.values()],
                'suffix' : None if self.sfxt is None else "compiled",
                'distmtx': self.dmtx}
