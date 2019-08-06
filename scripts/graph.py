import numpy as np
from Bio import Seq, SeqIO, SeqRecord
from block import Block

plus_strand = 1
minus_strand = -1

class Graph(object):
    """docstring for Graph"""
    def __init__(self):
        super(Graph, self).__init__()
        self.blocks = {}
        self.sequences = {}

    @classmethod
    def from_sequence(cls, sequence_name, sequence):
        new_graph = cls()
        b = Block.from_sequence(sequence_name, sequence)
        new_graph.blocks = {b.name:b}
        new_graph.sequences = {sequence_name: [(b.name, plus_strand)]}
        return new_graph

    @classmethod
    def fuse(cls, graph1, graph2):
        new_graph = Graph()
        new_graph.blocks = {}
        new_graph.blocks.update(graph1.blocks)
        new_graph.blocks.update(graph2.blocks)
        new_graph.sequences = {s:list(b) for s,b in list(graph1.sequences.items())+list(graph2.sequences.items())}
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
                for bi,b in enumerate(self.sequences[s]):
                    if b[0]==tmp_block.name:
                        orig_strand = b[1]
                        new_block_sequence = self.sequences[s][:bi] + [(x, orig_strand*y) for x,y in new_blocks] + self.sequences[s][bi+1:]
                        break
                self.sequences[s] = new_block_sequence

        remaining_blocks = set()
        for s in self.sequences:
            remaining_blocks.update([s[0] for s in self.sequences[s]])
        self.blocks = {b:self.blocks[b] for b in remaining_blocks}


    def extract(self, name):
        seq = ""
        for b,strand in self.sequences[name]:
            tmp_seq = self.blocks[b].extract(name)
            if strand==plus_strand:
                seq += tmp_seq
            else:
                seq += Seq.reverse_complement(tmp_seq)

        return seq


    def to_fasta(self, fname):
        SeqIO.write([SeqRecord.SeqRecord(seq=Seq.Seq("".join(c.consensus)), id=c.name, description='')
                    for c in self.blocks.values()], fname, format='fasta')


if __name__ == '__main__':
    g1 = Graph.from_sequence("S1", "xxxxxxACACACyyyy")
    g2 = Graph.from_sequence("S2", "zzzGTGTGTqqq")

    g = Graph.fuse(g1,g2)

    # this defines the homology hit format, meant to be parsed from paf
    hit = {'ref': {'name':g.sequences['S1'][0][0], 'start':6, 'end':12},
           'query': {'name':g.sequences['S2'][0][0], 'start':3, 'end':9},
           'cigar':'6M', 'orientation':minus_strand}

    g.merge_hit(hit)

    print(g.extract('S1'))
    print(g.extract('S2'))