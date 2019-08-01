import numpy as np
from cigar import Cigar

ref_complement = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}

def generate_block_name():
    alphabet = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    return "".join([alphabet[i] for i in np.random.randint(0,len(alphabet), 10)])

class Block(object):
    """docstring for Block"""
    def __init__(self):
        super(Block, self).__init__()
        self.sequences = {}
        self.consensus = None
        self.name = generate_block_name()

    @classmethod
    def from_sequence(cls, name, seq):
        new_block = cls()
        new_block.consensus = np.array(list(seq))
        new_block.sequences = {name:{}}
        return new_block


    @classmethod
    def from_alignment(cls, aln):
        new_block = cls()
        consensus = ""
        new_block.sequences = {}
        refpos = aln['ref_start']
        seqpos = 0
        R = {}
        S = {}
        tmpCigar = Cigar(aln["cigar"])
        for l,t in tmpCigar.items():
            consensus_pos = len(consensus)
            if t in ['S', 'H']:
                seqpos+=l
            elif t=='M':
                ref_block = np.array(list(aln['refseq'][refpos:refpos+l]))
                seq_block = np.array(list(aln['seq'][seqpos:seqpos+l]))
                diff = np.where(np.array(ref_block!=seq_block))[0]
                for i in diff:
                    S[i+consensus_pos] = seq_block[i]
                consensus += aln['refseq'][refpos:refpos+l]
                refpos += l
                seqpos += l
            elif t=='D':
                for i in range(l):
                    S[i+consensus_pos] = '-'
                consensus += aln['refseq'][refpos:refpos+l]
                refpos += l
            elif t=='I':
                for i in range(l):
                    R[i+consensus_pos] = '-'
                consensus += aln['seq'][seqpos:seqpos+l]
                seqpos += l

        new_block.sequences[aln['ref_name']] = R
        new_block.sequences[aln['seq_name']] = S
        new_block.consensus = np.array(list(consensus))

        return new_block

    @classmethod
    def from_cluster_alignment(cls, aln):
        new_block = cls()
        consensus = ""
        new_block.sequences = {}
        refpos = aln['ref_start']
        seqpos = 0
        consensus_pos = 0
        ref_map = [(refpos, 0)]
        seq_map = [(seqpos, 0)]

        # additional changes
        R = {}
        S = {}
        tmpCigar = Cigar(aln["cigar"])

        for l,t in tmpCigar.items():
            if t in ['S', 'H']:
                seqpos+=l
            elif t=='M':
                ref_block = np.array(list(aln['refseq'][refpos:refpos+l]))
                seq_block = np.array(list(aln['seq'][seqpos:seqpos+l]))
                diff = np.where(np.array(ref_block!=seq_block))[0]
                for i in diff:
                    S[i+consensus_pos] = seq_block[i]
                consensus+=aln['refseq'][refpos:refpos+l]
                refpos += l
                seqpos += l
            elif t=='D':
                for i in range(l):
                    S[i+consensus_pos] = '-'
                consensus+=aln['refseq'][refpos:refpos+l]
                refpos += l
            elif t=='I':
                for i in range(l):
                    R[i+consensus_pos] = '-'
                consensus+=aln['seq'][seqpos:seqpos+l]
                seqpos += l

            consensus_pos = len(consensus)

            ref_map.append((refpos, consensus_pos-refpos))
            seq_map.append((seqpos, consensus_pos-seqpos))

        ref_map = np.array(ref_map).T
        seq_map = np.array(seq_map).T

        for transform, extra_mods, sub_cluster in [(ref_map, R, aln['ref_cluster']),
                                                   (seq_map, S, aln['seq_cluster'])]:
            for s,muts in sub_cluster.items():
                old_pos = np.array(list(muts.keys()))
                new_pos = old_pos + transform[1][np.searchsorted(transform[0], old_pos, side='right')]
                tmp_seq = {newp:muts[oldp] for newp,oldp in zip(new_pos, old_pos)}
                tmp_seq.update(extra_mods)
                new_block.sequences[s] = tmp_seq

        new_block.consensus = np.array(list(consensus))

        return new_block

    def extract(self, seqname):
        tmp = np.copy(self.consensus)
        for p,s in self.sequences[seqname].items():
            tmp[p]=s

        return "".join(tmp[tmp!='-'])

    def copy(self):
        b = Block()
        b.consensus = np.copy(self.consensus)
        b.sequences = {s:{p:c for p.c in self.sequences[s].items()}
                       for s in self.sequences}
        return b

    def reverse_complement(self):
        from Bio import Seq
        new_block = Block()
        new_block.consensus = np.array(list(Seq.reverse_complement("".join(self.consensus))))
        n = len(self.consensus)
        for s in self.sequences:
            new_block.sequences[s] = {n-p: reverse_complement.get(c,c)
                                      for p,c in self.sequences[s].items()}
        return new_block


    def __getitem__(self, val):
        if isinstance(val, slice):
            b = Block()
            b.consensus = self.consensus[val.start:val.stop]
            for s, mods in self.sequences.items():
                b.sequences[s] = {p-val.start:c for p,c in self.sequences[s].items()
                                  if p>=val.start and p<val.stop}
            return b
        else:
            raise ValueError("item access not supported")

    def __len__(self):
        return len(self.consensus)


if __name__ == '__main__':
    seqs={'R1':"___ABCdEF123GHiJ", 'S1':'xxxABCD6789EFGHIJ',
          'R2':'ABCdEF123GHiJ', 'S2':'ABcdEFaaa23GHiJ'}
    aln = {'refseq':seqs["R1"], 'seq':seqs["S1"], 'cigar':'3S4M4I2M3D4M',
            'ref_name':'R1', 'seq_name':'S1', 'ref_start':3}

    b = Block.from_alignment(aln)

    aln = {'refseq':seqs["R2"], 'seq':seqs["S2"], 'cigar':'6M3I1D6M',
            'ref_name':'R2', 'seq_name':'S2', 'ref_start':0}

    c = Block.from_alignment(aln)

    caln = {"refseq":"".join(b.consensus), "seq":"".join(c.consensus),
            "ref_cluster":b.sequences, "seq_cluster":c.sequences,
            "ref_start":0, "cigar":"4M4D2M3I7M"}

    d  = Block.from_cluster_alignment(caln)

    for s in b.sequences:
        print(s, b.extract(s), seqs[s], b.extract(s)==seqs[s].strip('_x'))

    for s in c.sequences:
        print(s, c.extract(s), seqs[s], c.extract(s)==seqs[s].strip('_x'))


    for s in d.sequences:
        print(s, d.extract(s), seqs[s], d.extract(s)==seqs[s].strip('_x'))
