import os
import numpy as np

from utils import Strand, parsepaf
from block import Block
from Bio   import Seq, SeqIO, SeqRecord, Phylo

# ------------------------------------------------------------------------
# Global variables

outdir      = "data/graph"
maxselfmaps = 25

# ------------------------------------------------------------------------
# Helper functions

def asrecord(seq):
    return SeqRecord.SeqRecord(Seq.Seq(str(seq)), id=seq.name, name=seq.name)

def tryprint(msg, verbose):
    if verbose:
        print(msg)
    else:
        pass

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

    # --- Class methods ---

    @classmethod
    def fromseq(cls, name, seq):
        newg = cls()
        blk  = Block.fromseq(name, seq)
        newg.blks = {blk.id : blk}
        newg.seqs = {name : [(blk.id, Strand.Plus)]}
        newg.spos = {name :0}

        return newg

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
    def fromnwk(cls, path, seqs, save=True, verbose=False):
        print(path)
        T = Phylo.read(path, "newick")
        if T.count_terminals() == 1:
            return Graph()

        for n in T.get_terminals():
            seq      = asrecord(seqs[n.name])
            n.graph  = Graph.fromseq(seq.id, str(seq.seq).upper())
            n.name   = seq.id
            n.fapath = f"{Graph.blddir}/{n.name}.fasta"
            n.graph.tofasta()

        nnodes = 0
        for n in T.get_nonterminals(order='postorder'):
            # Simple graph "fuse". Straight concatenation
            nnodes += 1

            n.name   = f"NODE_{nnodes:07d}"
            n.graph  = Graph.fuse(n.clades[0].graph, n.clades[1].graph)
            n.fapath = os.path.join(*[Graph.blddir, f"{n.name}.fasta"])

            tryprint(f" -- node {n.name} with {n.count_terminals()} children", verbose)

            # Initial map from graph to itself
            n.graph, _ = n.graph.__mapandmerge(n.clades[0].fapath, n.clades[1].fapath,
                                    f"{Graph.blddir}/{n.name}")

            i, contin = 0, True
            while contin:
                tryprint(f"merge round {i}", verbose)
                n.graph.tofasta(f"{n.fapath}_iter_{i}")
                n.graph, contin = n.graph.__mapandmerge(f"{n.fapath}_iter_{i}", f"{n.fapath}_iter_{i}",
                                        f"{Graph.blddir}/{n.name}_iter_{i}")
                i += 1

                contin &= i < maxselfmaps

            tryprint(f"  --- Blocks: {len(n.graph.blks)}, length: {np.sum([len(b) for b in n.graph.blks.values()])}", verbose)

            n.graph.tofasta()

        G      = T.root.graph
        G.name = ".".join(os.path.basename(path).split(".")[:-1])

        if save:
            G.tojson()
            G.tofasta()

        return G

    # --- Internal methods ---

    def __mapandmerge(self, qpath, rpath, out):
        os.system(f"minimap2 -x asm5 -D -c {qpath} {rpath} 1>{out}.paf 2>log")

        paf = parsepaf(f"{out}.paf")
        paf.sort(key = lambda x: x['aligned_bases'], reverse=True)

        merged_blks = set()
        if len(paf) == 0:
            return self, False

        merged = False
        for hit in paf:
            if hit['qry']['name'] in merged_blks \
            or hit['ref']['name'] in merged_blks \
            or hit['ref']['name'] == hit['qry']['name'] \
            or hit['mapping_quality'] < 40:
                continue

            # citems = list(Cigar(hit['cigar']).items())

            self.merge(hit)
            merged = True
            merged_blks.add(hit['ref']['name'])
            merged_blks.add(hit['qry']['name'])

        self.rmempty()
        return self, merged

    # --- Instance methods ---

    def prune(self):
        blks_remain = set()
        for iso in self.seqs:
            blks_remain.update([s[0] for s in self.seqs[iso]])
        self.blks = {b:self.blks[b] for b in blks_remain}

        return

    def rmempty(self):
        for iso in self.seqs:
            goodblks = []
            popdblks = set()
            for i, (b, _) in enumerate(self.seqs[iso]):
                if b in popdblks:
                    continue

                bseq = self.blks[b].extract(iso)
                if bseq:
                    goodblks.append(i)
                else:
                    self.blks[b].muts.pop(iso)
                    popdblks.add(b)

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

        aln = {"ref_seq"     : "".join(refblk.seq),
               "qry_seq"     : "".join(qryblk.seq),
               "cigar"       : hit["cigar"],
               "ref_cluster" : refblk.muts,
               "qry_cluster" : qryblk.muts,
               "ref_start"   : hit["ref"]["start"],
               "qry_start"   : hit["qry"]["start"],
               "orientation" : hit["orientation"]}

        merged_blks, qryblks, refblks = Block.fromaln(aln)
        for merged_blk in merged_blks:
            self.blks[merged_blk.id] = merged_blk

        def update(blk, addblks, hit, strand):
            new_blks = []
            if hit['start'] > 0:
                left = blk[0:hit['start']]
                self.blks[left.id] = left
                new_blks.append((left.id, Strand.Plus))

            for b in addblks:
                new_blks.append((b.id, strand))

            if hit['end'] < len(blk):
                right = blk[hit['end']:]
                self.blks[right.id] = right
                new_blks.append((right.id, Strand.Plus))

            def replace(iso, newblks):
                new_blk_seq = []
                for b in self.seqs[iso]:
                    if b[0] == blk.id:
                        orig_strand = b[1]
                        if orig_strand == Strand.Plus:
                            tmp_new_blocks = [(x, orig_strand*y) for x,y in newblks]
                        else:
                            tmp_new_blocks = [(x, orig_strand*y) for x,y in newblks][::-1]
                        new_blk_seq.extend(tmp_new_blocks)
                    else:
                        new_blk_seq.append(b)

                self.seqs[iso] = new_blk_seq

            for iso in blk.muts.keys():
                replace(iso, new_blks)

        update(refblk_orig, refblks, hit['ref'], Strand.Plus)
        update(qryblk_orig, qryblks, hit['qry'], hit['orientation'])
        self.prune()

        return

    def extract(self, name, strip_gaps=True, verbose=False):
        seq = ""
        for (b, strand) in self.seqs[name]:
            tmp_seq = self.blks[b].extract(name, strip_gaps=False, verbose=verbose)
            if strand == Strand.Plus:
                seq += tmp_seq
            else:
                seq += Seq.reverse_complement(tmp_seq)

        start_pos = self.spos[name]

        if start_pos:
            seq = seq[start_pos:] + seq[:start_pos]

        if strip_gaps:
            seq = seq.replace('-', '')

        return seq

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
            aln = { J["Isolate_names"].index(s) :
                    self.blks[b].extract(s, strip_gaps=False) for s in self.blks[b].muts }
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
                edges[e] = {"Source":e[0], "Target":e[1], "Isolates":[pname]}

        for e in edges:
            nodes[e[0][0]]["Out_Edges"].append(edges[e])
            nodes[e[1][0]]["In_Edges"].append(edges[e])

        J["Nodes"] = nodes

        with open(f"{Graph.visdir}/{self.name}.json", 'w+') as fh:
            json.dump(J, fh)

        return

    def tofasta(self, path=None):
        if path is None:
            path = f"{Graph.alndir}/{self.name}.fasta"
        SeqIO.write([ SeqRecord.SeqRecord(seq=Seq.Seq("".join(c.seq)), id=c.id, description='')
                        for c in self.blks.values() ], path, format='fasta')

        return
