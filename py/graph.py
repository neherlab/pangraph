import os
import numpy as np

from glob  import glob
from utils import Strand, parsepaf, panic
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

def newstrand(s, t):
    if not isinstance(s, Strand) or not isinstance(t, Strand):
        raise TypeError(f"Expected an enum! Recieved {type(t)} and {type(s)}")

    if s != t:
        return Strand.Minus
    else:
        return Strand.Plus

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
        newg.name = name
        newg.blks = {blk.id : blk}
        newg.seqs = {name : [(blk.id, Strand.Plus, 0)]}
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

                            # import ipdb; ipdb.set_trace()
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
            return Graph()

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
            n.graph, _ = n.graph.__mapandmerge(n.clades[0].fapath, n.clades[1].fapath,
                                    f"{Graph.blddir}/{n.name}")

            i, contin = 0, True
            while contin:
                tryprint(f"----> merge round {i}", verbose)
                # check(seqs, T, n.graph)
                itrseq = f"{Graph.blddir}/{n.name}_iter_{i}"
                n.graph.tofasta(itrseq)
                n.graph, contin = n.graph.__mapandmerge(itrseq, itrseq, f"{Graph.blddir}/{n.name}_iter_{i}")
                i += 1

                contin &= i < maxselfmaps

            tryprint(f"-- Blocks: {len(n.graph.blks)}, length: {np.sum([len(b) for b in n.graph.blks.values()])}\n", verbose)

            n.graph.tofasta(f"{Graph.blddir}/{n.name}")

        G      = T.root.graph
        G.name = ".".join(os.path.basename(path).split(".")[:-1])
        check(seqs, T, G)

        if save:
            G.tojson()
            G.tofasta()

        return G

    @classmethod
    def cleanbld(cls):
        for p in glob(f"{cls.blddir}/*"):
            os.remove(p)

    # --- Internal methods ---

    def __mapandmerge(self, qpath, rpath, out):
        os.system(f"minimap2 -x asm5 -D -c {qpath}.fasta {rpath}.fasta 1>{out}.paf 2>log")

        paf = parsepaf(f"{out}.paf")
        paf.sort(key = lambda x: x['aligned_bases'], reverse=True)

        merged_blks = set()
        if len(paf) == 0:
            return self, False

        merged = False
        for hit in paf:
            if hit['qry']['name'] in merged_blks \
            or hit['ref']['name'] in merged_blks \
            or(hit['ref']['name'] <= hit['qry']['name'] and qpath == rpath) \
            or hit['mapping_quality'] < 40:
                continue

            # TODO: We originally had an if statement here that checked for mapping quality.
            # This is gone, although replaced with the cutoff above.
            # Think through if we would like to do something different here.

            # citems = list(Cigar(hit['cigar']).items())

            tryprint(f"------> merge {hit['ref']['name']} with {hit['qry']['name']}", False)
            self.merge(hit)
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

        # if "FYLWYFIIUV" in self.blks:
        #     import ipdb; ipdb.set_trace()

        return

    def purgeempty(self):
        for iso in self.seqs:
            goodblks = []
            popblks  = set()
            for i, (b, _, n) in enumerate(self.seqs[iso]):
                # if b == "BDWFSBDNCP":
                #     import ipdb; ipdb.set_trace()

                if b in popblks:
                    continue

                bseq = self.blks[b].extract(iso, n)
                if bseq:
                    goodblks.append(i)
                else:
                    self.blks[b].muts.pop((iso, n))

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

        aln = {"ref_seq"     : "".join(refblk.seq),
               "qry_seq"     : "".join(qryblk.seq),
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
            # if blk.id == "TCSGYOLBUN":
            #     pos   = [0]
            #     blks  = self.seqs['carb003_3']
            #     for b, strand, num in blks:
            #         pos.append(pos[-1] + len(self.blks[b].extract('carb003_3', num)))
            #     pos = pos[1:]

            #     # import ipdb; ipdb.set_trace()

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

            # for nb in new_blks:
            #     if nb[0] == 'IIAAMUUZVS':
            #         import ipdb; ipdb.set_trace()
            #         iso = "carb003_3"
            #         pos = [0]
            #         for b in self.seqs[iso]:
            #             pos.append(pos[-1] + len(self.blks[b[0]].extract(iso, b[2])))
            #         pos = pos[1:]

            def replace(tag, blk, newblks):
                iso = tag[0]

                new_blk_seq = []
                for b in self.seqs[iso]:
                    if b[0] == blk.id and b[2] == tag[1]:
                        # if iso == "carb003_3" and b[0] == "TCSGYOLBUN" and tag[1] == 2:
                        #     import ipdb; ipdb.set_trace()

                        orig_strand    = b[1]
                        tmp_new_blocks = [(ID, newstrand(orig_strand, ns), isomap[blk.id][tag][1]) if merged else (ID, newstrand(orig_strand, ns), b[2]) for ID, ns, merged in newblks]
                        if orig_strand == Strand.Minus:
                            tmp_new_blocks = tmp_new_blocks[::-1]

                        new_blk_seq.extend(tmp_new_blocks)
                    else:
                        new_blk_seq.append(b)

                # if iso == "carb003_3" and blk.id == "TCSGYOLBUN" and tag[1] == 2:
                #     import ipdb; ipdb.set_trace()
                self.seqs[iso] = new_blk_seq

            for tag in blk.muts.keys():
                replace(tag, blk, new_blks)

        # print(f"UPDATING {refblk_orig.id} and {qryblk_orig.id}")
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

        print(f"----> Saving to {Graph.visdir}/{self.name}.json")
        with open(f"{Graph.visdir}/{self.name}.json", 'w+') as fh:
            json.dump(J, fh)

        return

    def tofasta(self, path=None):
        if path is None:
            path = f"{Graph.alndir}/{self.name}.fasta"
        elif not path.endswith(".fasta"):
            path = path + ".fasta"
        SeqIO.write(sorted([ SeqRecord.SeqRecord(seq=Seq.Seq("".join(c.seq)), id=c.id, description='')
            for c in self.blks.values() ], key=lambda x: len(x), reverse=True), path, format='fasta')

        return
