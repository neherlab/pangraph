import numpy as np
from Bio import Seq
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
        orig_ref_block = self.blocks[hit[0][0]]
        orig_seq_block = self.blocks[hit[1][0]]
        ref_block = orig_ref_block[hit[0][1]:hit[0][2]]
        seq_block = orig_seq_block[hit[1][1]:hit[1][2]]
        if hit[2]["orientation"]==minus_strand:
            seq_block = seq_block.reverse_complement()

        aln = {"refseq":"".join(ref_block.consensus),
               "seq":"".join(seq_block.consensus),
               "cigar":hit[2]["cigar"], "ref_cluster":ref_block.sequences,
               "seq_cluster":seq_block.sequences, 'ref_start':0}

        merged_block = Block.from_cluster_alignment(aln)
        self.blocks[merged_block.name] = merged_block

        for tmp_block, subhit, strand in [(orig_ref_block, hit[0], plus_strand),
                                          (orig_seq_block, hit[1], hit[2]["orientation"])]:
            new_blocks = []
            if subhit[1]:
                left = tmp_block[0:subhit[1]]
                new_blocks.append((left.name, plus_strand))
                self.blocks[left.name] = left
            new_blocks.append((merged_block.name, strand))
            if subhit[2]<len(tmp_block):
                right = tmp_block[subhit[2]:]
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

        ## replace block in original sequences by new blocks
        # ideally do recursively in that every block can be a reconstructed from its constituents


    def extract(self, name):
        seq = ""
        for b,strand in self.sequences[name]:
            tmp_seq = self.blocks[b].extract(name)
            print(b,strand,tmp_seq)
            if strand==plus_strand:
                seq += tmp_seq
            else:
                seq += Seq.reverse_complement(tmp_seq)

        return seq

if __name__ == '__main__':
    g1 = Graph.from_sequence("S1", "xxxxxxACACACyyyy")
    g2 = Graph.from_sequence("S2", "zzzGTGTGTqqq")

    g = Graph.fuse(g1,g2)

    hit = [(g.sequences['S1'][0][0], 6, 12),
           (g.sequences['S2'][0][0], 3, 9),
           {'cigar':'6M', 'orientation':minus_strand}]

    g.merge_hit(hit)
