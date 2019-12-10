from collections import defaultdict
import numpy as np

from Bio   import Seq, SeqIO, SeqRecord
from block import Block

# ------------------------------------------------------------------------
# Globals

plus_strand  = +1
minus_strand = -1

def panic(msg):
    raise ValueError(f"Panic: {msg}")

# ------------------------------------------------------------------------
# Graph class

class Graph(object):
    """docstring for Graph"""
    def __init__(self):
        super(Graph, self).__init__()
        self.blks = {}
        self.seqs = {}
        self.spos = {}
        self.ST = None

    @classmethod
    def from_seq(cls, sequence_name, sequence):
        newg = cls()
        blk  = Block.from_seq(sequence_name, sequence)
        newg.blks = {blk.id:blk}
        newg.seqs = {sequence_name: [(blk.id, plus_strand)]}
        newg.spos = {sequence_name:0}

        return newg

    @classmethod
    def fuse(cls, graph1, graph2):
        newg = Graph()
        newg.blks = {}
        newg.blks.update(graph1.blks)
        newg.blks.update(graph2.blks)
        newg.seqs = {s:list(b) for s,b in list(graph1.seqs.items())+list(graph2.seqs.items())}
        newg.spos = {s:b for s,b in list(graph1.spos.items())+list(graph2.spos.items())}

        return newg

    def merge_hit(self, hit):
        # print("-------> starting to merge hit")
        refblk_orig = self.blks[hit['ref']['name']]
        qryblk_orig = self.blks[hit['qry']['name']]
        # print("-------> grabbed blocks")

        # As we slice here, we DONT need to remember the starting position.
        # This is why in from_aln(aln) we set the start index to 0
        refblk = refblk_orig[hit['ref']['start']:hit['ref']['end']]
        qryblk = qryblk_orig[hit['qry']['start']:hit['qry']['end']]

        if hit["orientation"] == minus_strand:
            qryblk = qryblk.rev_cmpl()

        aln = {"ref_seq"     : "".join(refblk.seq),
               "qry_seq"     : "".join(qryblk.seq),
               "cigar"       : hit["cigar"],
               "ref_cluster" : refblk.muts,
               "qry_cluster" : qryblk.muts,
               "ref_start"   : hit["ref"]["start"],
               "qry_start"   : hit["qry"]["start"],
               "orientation" : hit["orientation"]}

        # print("CIGAR:", hit["cigar"])
        # import ipdb
        # ipdb.set_trace()

        # print("-------> building merged block")
        merged_blks, qryblks, refblks = Block.from_aln(aln)
        for merged_blk in merged_blks:
            self.blks[merged_blk.id] = merged_blk

        # if len(merged_blks) == 1 and len(qryblks) == 1 and len(refblks) == 1:
            # import ipdb
            # ipdb.set_trace()

        def update(blk, addblks, hit, strand):
            new_blks = []
            if hit['start'] > 0:
                left = blk[0:hit['start']]
                self.blks[left.id] = left
                new_blks.append((left.id, plus_strand))

            for b in addblks:
                new_blks.append((b.id, strand))

            if hit['end'] < len(blk):
                right = blk[hit['end']:]
                self.blks[right.id] = right
                new_blks.append((right.id, plus_strand))

            def replace(iso, newblks):
                new_blk_seq = []
                # last_blk    = 0
                for bi, b in enumerate(self.seqs[iso]):
                    if b[0] == blk.id:
                        orig_strand = b[1]
                        if orig_strand == plus_strand:
                            tmp_new_blocks = [(x, orig_strand*y) for x,y in newblks]
                        else:
                            tmp_new_blocks = [(x, orig_strand*y) for x,y in newblks][::-1]
                        new_blk_seq.extend(tmp_new_blocks)
                        # NOTE: This break assumes there is one block instance per genome
                        # new_blk_seq = self.seqs[iso][last_blk:bi] + tmp_new_blocks
                        # last_blk    = bi+1
                        # break 
                    else:
                        new_blk_seq.append(b)
                # new_blk_seq += self.seqs[iso][last_blk:]

                # if len(new_blk_seq) > 0:
                self.seqs[iso] = new_blk_seq

            # print(f"------> replacing block '{blk.id}' with new blocks")
            for iso in blk.muts.keys():
                replace(iso, new_blks)

        update(refblk_orig, refblks, hit['ref'], plus_strand)
        update(qryblk_orig, qryblks, hit['qry'], hit['orientation'])
        # print("------> pruning bad blocks")
        self.prune_blocks()

        return

    def prune_blocks(self):
        blks_remain = set()
        # ---- Remove when done debugging ----- #
        old_blks = set([ b[0] for s in self.seqs for b in self.seqs[s] ])
        # ------------------------------------- #
        for s in self.seqs:
            blks_remain.update([s[0] for s in self.seqs[s]])
        # ---- Remove when done debugging ----- #
        # print(f"------> Removed blocks {old_blks.difference(blks_remain)}")
        # ------------------------------------- #
        self.blks = {b:self.blks[b] for b in blks_remain}

        return

    def extract(self, name, strip_gaps=True, verbose=False):
        seq = ""
        for i, (b, strand) in enumerate(self.seqs[name]):
            # print("***", name, i, b, strand, sum(1 for n in seq if n != '-'))
            tmp_seq = self.blks[b].extract(name, strip_gaps=False, verbose=verbose)
            # print("----->", "".join(tmp_seq)[0:100], "...")
            if strand==plus_strand:
                seq += tmp_seq
            else:
                seq += Seq.reverse_complement(tmp_seq)
        start_pos = self.spos[name]
        # print(name, start_pos)
        if start_pos:
            seq = seq[start_pos:] + seq[:start_pos]

        if strip_gaps:
            seq = seq.replace('-', '')

        return seq


    def prune_empty(self):
        for seq in self.seqs:
            good_blocks = []
            popped_blks = set()
            for bi, (b, _) in enumerate(self.seqs[seq]):
                if b in popped_blks:
                    continue

                bseq = self.blks[b].extract(seq)
                if bseq:
                    good_blocks.append(bi)
                else:
                    print(f"Popping block {b} from {seq}. Obtained bseq")
                    self.blks[b].muts.pop(seq)
                    popped_blks.add(b)
            self.seqs[seq] = [self.seqs[seq][bi] for bi in good_blocks]

        return


    def flip_edge(self, edge):
        return ((edge[1][0], minus_strand*edge[1][1]), (edge[0][0], minus_strand*edge[0][1]))

    def get_edges(self):
        edges = defaultdict(list)
        for seq, p in self.seqs.items():
            for bi,b in enumerate(p):
                if p[bi-1][0]<b[0]:
                    label = (p[bi-1], b)
                else:
                    label = self.flip_edge((p[bi-1], b))

                edges[label].append(seq)
        return edges


    def prune_transitive_edges(self):
        edges = self.get_edges()
        transitive = []
        # to label an edge, we need pairs of blocks together with their orientation
        # edges have an inversion symmetry: flipping order and strand results in the same edge
        # this is already standardized by the self.get_edges() such that each edge starts with
        # the alphabetically earlier block
        for (b1,s1),(b2,s2) in edges:
            if self.blks[b1].muts.keys()==self.blks[b2].muts.keys() \
                and set(edges[((b1,s1),(b2,s2))])==self.blks[b1].muts.keys():
                transitive.append([(b1,s1),(b2,s2)])

        chains = {}
        # bs1, bs2 are (block, strand) pairs
        # the chaining needs
        for (b1, s1), (b2, s2) in transitive:
            if b1 in chains and b2 in chains:
                c1 = chains[b1]
                c2 = chains[b2]
                if c1==c2: # this would circularize/duplicate the chain
                    continue
                if (b1,s1)==c1[-1] and (b2,s2)==c2[0]:
                    new_chain = c1 + c2
                elif (b1,s1)==c1[-1] and (b2,s2*minus_strand)==c2[-1]:
                    new_chain = c1 + [(b,s*minus_strand) for b,s in c2[::-1]]
                elif (b1,s1*minus_strand)==c1[0] and (b2,s2*minus_strand)==c2[-1]:
                    new_chain = c2 + c1
                elif (b1,s1*minus_strand)==c1[0] and (b2,s2)==c2[0]:
                    new_chain = [(b,s*minus_strand) for b,s in c1[::-1]] + c2
                else:
                    print("not covered:", b1, s1, b2, s1, c1, c2)
                for b,s in new_chain:
                    chains[b] = new_chain
            elif b1 in chains:
                if (b1,s1)==chains[b1][-1]:
                    chains[b1].append((b2,s2))
                elif (b1, s1*minus_strand)==chains[b1][0]:
                    chains[b1].insert(0, (b2,s2*minus_strand))
                else:
                    panic("chains should be unbranched")
                chains[b2] = chains[b1]
            elif b2 in chains:
                if (b2,s2*minus_strand)==chains[b2][-1]:
                    chains[b2].append((b1, s1*minus_strand))
                elif (b2,s2)==chains[b2][0]:
                    chains[b2].insert(0,(b1,s1))
                else:
                    panic("chains should be unbranched")
                chains[b1] = chains[b2]
            else:
                chains[b1] = [(b1,s1),(b2,s2)]
                chains[b2] = chains[b1]

        unique_chains = {id(x):x for x in chains.values()}.values()
        for c in unique_chains:
            seqs = set.union(set(edges[(c[0],c[1])]), set(edges[self.flip_edge((c[0],c[1]))]))

            c_start, c_end = c[0], c[-1]
            concat = Block.cat([self.blks[b[0]] \
                        if b[1]==plus_strand \
                        else self.blks[b[0]].rev_cmpl() for b in c])
            for seq in seqs:
                # - find first/last block
                # - use strand in seq and chain to determine orientation
                # - fwd and end>begin -> internal
                # - rev  and begin<end -> internal
                # - otherwise -> wraps around the origin

                p = self.seqs[seq]
                p_blocks = [b[0] for b in p] #sequence regardless of strand
                ci_start = p_blocks.index(c_start[0])
                ci_end = p_blocks.index(c_end[0])

                strand = plus_strand if p[ci_start]==c_start else minus_strand
                p_start = ci_start if strand==plus_strand else ci_end
                p_end = ci_end if strand==plus_strand else ci_start
                if p_start<p_end: # internal segment
                    self.seqs[seq] = p[:p_start] + [(concat.id, strand)] + p[p_end+1:]
                else:
                    self.seqs[seq] = [(concat.id, strand)] + p[p_end+1:p_start]
                    self.spos[seq] = self.spos.get(seq, 0) + np.sum([len(self.blks[b[0]]) for b in p[p_start:]])

            self.blks[concat.id] = concat
            # print(self.blocks.keys())
            for b in c:
                self.blks.pop(b[0])


    def shared_block_length(self):
        from itertools import combinations
        bl = defaultdict(list)
        for b in self.blks.values():
            l = len(b)
            for s1, s2 in combinations(b.muts, r=2):
                bl[(s1,s2)].append(l)
                bl[(s2,s1)].append(l)

        return bl


    def sub_graph(self, seqs):
        newg = Graph()
        newg.seqs = {s:list(self.seqs[s]) for s in seqs}
        newg.spos = {s:self.spos[s] for s in seqs}
        newg.blks = {}
        blks_remain = set()
        for s in newg.seqs:
            blks_remain.update([x[0] for x in newg.seqs[s]])

        for bname in blks_remain:
            newblk = self.blks[bname].copy(keep_name=True)
            newblk.muts = {s:v for s,v in newblk.muts.items() if s in seqs}
            newg.blks[bname] = newblk

        newg.prune_blocks()
        return newg


    def make_suffix_tree(self):
        # make a generalized suffix tree containing 5-3 and 3-5 of all sequences
        # terminate at sequence_name_fwd/rev_pos
        #
        from suffix_tree import Tree
        strings = {s+'_fwd':seq+seq[:-1] for s, seq in self.seqs.items()}
        strings.update({s+'_rev':[(b, s*minus_strand) for b,s in seq[::-1]+seq[::-1][:-1]]
                                   for s, seq in self.seqs.items()})
        # print(strings)
        self.ST = Tree(strings)


    def strip_suffix_tree(self):
        """
        we added the sequences twice, hence we only retain suffixes that are at most the
        length of the longest sequence in a particular branch of the tree.
        this doesn't actually work -- one needs to prune more specifically suffixes that don't
        start in copy one
        """
        seq_length = {s:len(seq) for s, seq in self.seqs.items()}
        def assign_max(x):
            if x.is_leaf():
                x.max_length = seq_length[x.str_id[:-4]]
            else:
                x.max_length = np.max([yv.max_length for yv in x.children.values()])
        self.ST.root.post_order(assign_max)

        def prune_redundant(x):
            if x.is_internal():
                x.children = {k:v for k,v in x.children.items() if v.string_depth()<=v.max_length+2}
        self.ST.root.pre_order(prune_redundant)


    def find_common_substrings(self, seqs):
        common_substrings = {}
        seq_len = {s:len(self.seqs[s]) for s in self.seqs}
        def attach_tip_names(x):
            if x.is_leaf():
                x.tip_names = {x.str_id[:-4]}
                x.tip_names_pos = {(x.str_id, (x.path.end - x.string_depth()) % seq_len[x.str_id[:-4]])}
            else:
                x.tip_names = set()
                x.tip_names_pos = set()
                for c in x.children.values():
                    x.tip_names_pos.update(c.tip_names_pos)
                    x.tip_names.update(c.tip_names)
        self.ST.root.post_order(attach_tip_names)

        def get_common_substrings(x):
            if len(x.path) and len(seqs.intersection(x.tip_names)) == len(seqs):
                label = tuple(sorted([s for s in x.tip_names_pos if s[0][:-4] in seqs]))
                if label in common_substrings:
                    if len(common_substrings[label])<len(x.path):
                        common_substrings[label] = x.path
                else:
                    common_substrings[label] = x.path

        self.ST.root.pre_order(get_common_substrings)
        covered = { s : np.zeros(len(self.seqs[s]), dtype=int) for s in seqs    }
        for cs in common_substrings:
            l = len(common_substrings[cs])
            for s,p in cs:
                lseq = seq_len[s[:-4]]
                if s.endswith('fwd'):
                    start = p
                    end = start + l
                else:
                    end = (2*lseq-1) - p + 1
                    start = end - l

                start_mod = start%lseq
                end_mod = end%lseq or lseq
                # print(l, lseq, p, start, end, start_mod, end_mod)
                if start_mod<end_mod:
                    covered[s[:-4]][start_mod:end_mod] += 1
                else:
                    covered[s[:-4]][start_mod:] += 1
                    covered[s[:-4]][:end_mod] += 1

        return common_substrings, covered

    def to_fasta(self, fname):
        SeqIO.write([SeqRecord.SeqRecord(seq=Seq.Seq("".join(c.seq)), id=c.id, description='')
                    for c in self.blks.values()], fname, format='fasta')


    def to_json(self, fname, min_length=500):
        J = {}
        cleaned_seqs = {s:[b for b in self.seqs[s] if len(self.blks[b[0]])>min_length]
                             for s in self.seqs}
        relevant_blocks = set()
        for s in cleaned_seqs.values():
            relevant_blocks.update([b[0] for b in s])
        # import ipdb; ipdb.set_trace()

        J['Isolate_names'] = list(cleaned_seqs.keys())
        J['Plasmids'] = [[x for x in cleaned_seqs[s]] for s in J['Isolate_names']]
        nodes = {}
        for b in relevant_blocks:
            nodes[b] = {"ID":b,
                        "Genomes": {"Consensus":''.join(self.blks[b].seq),
                        "Alignment": {J["Isolate_names"].index(s):self.blks[b].extract(s, strip_gaps=False)
                                                                 for s in self.blks[b].muts}},
                        "Out_Edges":[], "In_Edges":[]}

        edges = {}
        for pname, p in zip(range(len(J["Isolate_names"])), J['Plasmids']):
            for i in range(len(p)-1):
                e = (p[i], p[i+1])
                if e in edges:
                    edges[e]["Isolates"].append(pname)
                else:
                    edges[e] = {"Source":e[0], "Target":e[1], "Isolates":[pname]}
            e = (p[-1], p[0])
            if e in edges:
                edges[e]["Isolates"].append(pname)
            else:
                edges[e] = {"Source":e[0], "Target":e[1], "Isolates":[pname]}
        for e in edges:
            nodes[e[0][0]]["Out_Edges"].append(edges[e])
            nodes[e[1][0]]["In_Edges"].append(edges[e])
        J["Nodes"] = nodes
        import json
        with open(fname, 'w') as fh:
            json.dump(J, fh)
