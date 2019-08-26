import numpy as np
from Bio import Seq, SeqIO, SeqRecord
from collections import defaultdict
from block import Block

plus_strand = 1
minus_strand = -1

class Graph(object):
    """docstring for Graph"""
    def __init__(self):
        super(Graph, self).__init__()
        self.blocks = {}
        self.sequences = {}
        self.sequence_start = {}

    @classmethod
    def from_sequence(cls, sequence_name, sequence):
        new_graph = cls()
        b = Block.from_sequence(sequence_name, sequence)
        new_graph.blocks = {b.name:b}
        new_graph.sequences = {sequence_name: [(b.name, plus_strand)]}
        new_graph.sequence_start = {sequence_name:0}
        return new_graph

    @classmethod
    def fuse(cls, graph1, graph2):
        new_graph = Graph()
        new_graph.blocks = {}
        new_graph.blocks.update(graph1.blocks)
        new_graph.blocks.update(graph2.blocks)
        new_graph.sequences = {s:list(b) for s,b in list(graph1.sequences.items())+list(graph2.sequences.items())}
        new_graph.sequence_start = {s:b for s,b in list(graph1.sequence_start.items())+list(graph2.sequence_start.items())}
        return new_graph

    def merge_hit(self,hit):
        orig_ref_block = self.blocks[hit['ref']['name']]
        orig_seq_block = self.blocks[hit['query']['name']]
        ref_block = orig_ref_block[hit['ref']['start']:hit['ref']['end']]
        seq_block = orig_seq_block[hit['query']['start']:hit['query']['end']]

        if hit["orientation"]==minus_strand:
            seq_block = seq_block.reverse_complement()

        aln = {"ref_seq":"".join(ref_block.consensus),
               "query_seq":"".join(seq_block.consensus),
               "cigar":hit["cigar"], "ref_cluster":ref_block.sequences,
               "query_cluster":seq_block.sequences, 'ref_start':0}

        merged_block = Block.from_cluster_alignment(aln)
        self.blocks[merged_block.name] = merged_block

        for tmp_block, subhit, strand in [(orig_ref_block, hit['ref'], plus_strand),
                                          (orig_seq_block, hit['query'], hit["orientation"])]:
            new_blocks = []
            if subhit['start']:
                left = tmp_block[0:subhit['start']]
                new_blocks.append((left.name, plus_strand))
                self.blocks[left.name] = left
            new_blocks.append((merged_block.name, strand))
            if subhit['end']<len(tmp_block):
                right = tmp_block[subhit['end']:]
                new_blocks.append((right.name, plus_strand))
                self.blocks[right.name] = right

            for s in tmp_block.sequences:
                new_block_sequence = []
                last_block = 0
                for bi,b in enumerate(self.sequences[s]):
                    if b[0]==tmp_block.name:
                        orig_strand = b[1]
                        if orig_strand==plus_strand:
                            tmp_new_blocks = [(x, orig_strand*y) for x,y in new_blocks]
                        else:
                            tmp_new_blocks = [(x, orig_strand*y) for x,y in new_blocks][::-1]
                        new_block_sequence = self.sequences[s][last_block:bi] + tmp_new_blocks
                        last_block = bi+1
                        # print(b[0], bi, len(new_block_sequence), len(self.sequences[s]))
                        break
                new_block_sequence += self.sequences[s][last_block:]
                if len(new_block_sequence):
                    self.sequences[s] = new_block_sequence

        self.prune_blocks()

    def prune_blocks(self):
        remaining_blocks = set()
        for s in self.sequences:
            remaining_blocks.update([s[0] for s in self.sequences[s]])
        self.blocks = {b:self.blocks[b] for b in remaining_blocks}


    def extract(self, name, strip_gaps=True):
        seq = ""
        for b,strand in self.sequences[name]:
            tmp_seq = self.blocks[b].extract(name, strip_gaps=False)
            if strand==plus_strand:
                seq += tmp_seq
            else:
                seq += Seq.reverse_complement(tmp_seq)
        start_pos = self.sequence_start[name]
        if start_pos:
            # import ipdb; ipdb.set_trace()
            seq = seq[start_pos:] + seq[:start_pos]

        if strip_gaps:
            seq = seq.replace('-', '')

        return seq

    def prune_empty(self):
        for seq in self.sequences:
            good_blocks = []
            for bi, (b,s) in enumerate(self.sequences[seq]):
                bseq = self.blocks[b].extract(seq)
                if bseq:
                    good_blocks.append(bi)
                else:
                    self.blocks[b].sequences.pop(seq)
                    print("Pop", seq, b)
            self.sequences[seq] = [self.sequences[seq][bi] for bi in good_blocks]

    def get_edges(self):
        edges = defaultdict(list)
        for seq, p in self.sequences.items():
            for bi,b in enumerate(p):
                if p[bi-1][0]<b[0]:
                    label = (p[bi-1], b)
                else:
                    label = ((b[0], minus_strand*b[1]), (p[bi-1][0], minus_strand*p[bi-1][1]))

                edges[label].append(seq)
        return edges


    def prune_transitive_edges(self):
        edges = self.get_edges()
        transitive = []
        # to label an edge, we need pairs of blocks together with their orientation
        # edges have an inversion symmetry: flipping order and strand results in the same edge
        # this is already standardized by the self.get_edges() such that each edge starts with
        # the alphabetically earlier block
        #
        # The latter needs to check that there are no intervening edges.
        for (b1,s1),(b2,s2) in edges:
            if self.blocks[b1].sequences.keys()==self.blocks[b2].sequences.keys() \
                and set(edges[((b1,s1),(b2,s2))])==self.blocks[b1].sequences.keys():
                transitive.append([(b1,s1),(b2,s2)])

        chains = {}
        # bs1, bs2 are (block, strand) pairs
        for bs1,bs2 in transitive:
            if bs1 in chains:
                if bs1==chains[bs1][-1]:
                    chains[bs1].append(bs2)
                elif bs1==chains[bs1][0]:
                    chains[bs1].insert(0,bs2)
                else:
                    raise ValueError("chains should be unbranched")
                chains[bs2] = chains[bs1]
            elif b2 in chains:
                if bs2==chains[bs2][-1]:
                    chains[b2].append(b1)
                elif bs2==chains[bs2][0]:
                    chains[bs2].insert(0,bs1)
                else:
                    raise ValueError("chains should be unbranched")
                chains[bs1] = chains[bs2]
            else:
                chains[bs1] = [bs1,bs2]
                chains[bs2] = chains[bs1]

        unique_chains = {id(x):x for x in chains.values()}.values()
        print("Chains:", unique_chains)
        for c in unique_chains:
            c_start,c_end = c[0],c[-1]
            if len(c)>3:
                import ipdb; ipdb.set_trace()

            concat = Block.concatenate([self.blocks[b[0]] if b[1]==plus_strand
                                        else self.blocks[b[0]].reverse_complement()
                                        for b in c])
            for seq in edges[(c[0],c[1])]:
                # alternative
                # - find first/last block
                # - use strand in seq and chain to determine orientation
                # - fwd and end>begin -> internal
                # - rev  and begin<end -> internal
                # - otherwise -> wraps around the origin
                #
                p = self.sequences[seq]
                p_blocks = [b[0] for b in p] #sequence regardless of strand
                ci_start = p_blocks.index(c_start[0])
                ci_end = p_blocks.index(c_end[0])

                strand = plus_strand if p[ci_start]==c_start else minus_strand
                p_start = ci_start if strand==plus_strand else ci_end
                p_end = ci_end if strand==plus_strand else ci_start
                if p_start<p_end: # internal segment
                    self.sequences[seq] = p[:p_start] + [(concat.name, strand)] + p[p_end+1:]
                else:
                    self.sequences[seq] = [(concat.name, strand)] + p[p_end+1:p_start]
                    self.sequence_start[seq] = np.sum([len(self.blocks[b[0]]) for b in p[p_start:]])

            self.blocks[concat.name] = concat
            for b in c:
                self.blocks.pop(b[0])

        return transitive


    def to_fasta(self, fname):
        SeqIO.write([SeqRecord.SeqRecord(seq=Seq.Seq("".join(c.consensus)), id=c.name, description='')
                    for c in self.blocks.values()], fname, format='fasta')


    def to_json(self, fname, min_length=500):
        J = {}
        cleaned_sequences = {s:[b for b in self.sequences[s] if len(self.blocks[b[0]])>min_length]
                             for s in self.sequences}
        relevant_blocks = set()
        for s in cleaned_sequences.values():
            relevant_blocks.update([b[0] for b in s])
        # import ipdb; ipdb.set_trace()

        J['Isolate_names'] = list(cleaned_sequences.keys())
        J['Plasmids'] = [[x[0] for x in cleaned_sequences[s]] for s in J['Isolate_names']]
        nodes = {}
        for b in relevant_blocks:
            nodes[b] = {"ID":b,
                        "Genomes":{"Consensus":''.join(self.blocks[b].consensus),
                        "Alignment":{J["Isolate_names"].index(s):self.blocks[b].extract(s, strip_gaps=False)
                                                                 for s in self.blocks[b].sequences}},
                        "Out_Edges":[], "In_Edges":[]}

        edges = {}
        for pname, p in zip(range(len(J["Isolate_names"])), J['Plasmids']):
            for i in range(len(p)-1):
                e = (p[i],p[i+1])
                if e in edges:
                    edges[e]["Isolates"].append(pname)
                else:
                    edges[e] = {"Source":e[0], "Target":e[1], "Isolates":[pname]}
            e = (p[-1],p[0])
            if e in edges:
                edges[e]["Isolates"].append(pname)
            else:
                edges[e] = {"Source":e[0], "Target":e[1], "Isolates":[pname]}
        for e in edges:
            nodes[e[0]]["Out_Edges"].append(edges[e])
            nodes[e[1]]["In_Edges"].append(edges[e])
        J["Nodes"] = nodes
        import json
        with open(fname, 'w') as fh:
            json.dump(J, fh)


    def sub_graph(self, sequences):
        new_graph = Graph()
        new_graph.sequences = {s:list(self.sequences[s]) for s in sequences}
        new_graph.sequence_start = {s:self.sequence_start[s] for s in sequences}
        new_graph.blocks = {}
        for bname,b in self.blocks.items():
            new_block = b.copy(keep_name=True)
            new_block.sequences = {s:v for s,v in new_block.sequences.items()
                                                  if s in sequences}
            new_graph.blocks[bname] = new_block

        new_graph.prune_blocks()
        return new_graph



if __name__ == '__main__':
    g1 = Graph.from_sequence("S1", "xxxxxxACACACyyyy")
    g2 = Graph.from_sequence("S2", "zzzGTCTGTqqq")

    g = Graph.fuse(g1,g2)

    # this defines the homology hit format, meant to be parsed from paf
    hit = {'ref': {'name':g.sequences['S1'][0][0], 'start':6, 'end':12},
           'query': {'name':g.sequences['S2'][0][0], 'start':3, 'end':9},
           'cigar':'6M', 'orientation':minus_strand}

    g.merge_hit(hit)

    print(g.extract('S1'))
    print(g.extract('S2'))
