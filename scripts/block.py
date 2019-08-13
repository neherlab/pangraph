import numpy as np
from cigar import Cigar

reverse_complement = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}

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
        ref_pos = aln['ref_start']
        query_pos = 0
        R = {}
        S = {}
        tmpCigar = Cigar(aln["cigar"])
        for l,t in tmpCigar.items():
            consensus_pos = len(consensus)
            if t in ['S', 'H']:
                query_pos+=l
            elif t=='M':
                ref_block = np.array(list(aln['ref_seq'][ref_pos:ref_pos+l]))
                query_block = np.array(list(aln['query_seq'][query_pos:query_pos+l]))
                diff = np.where(np.array(ref_block!=query_block))[0]
                for i in diff:
                    S[i+consensus_pos] = query_block[i]
                consensus += aln['ref_seq'][ref_pos:ref_pos+l]
                ref_pos += l
                query_pos += l
            elif t=='D':
                for i in range(l):
                    S[i+consensus_pos] = '-'
                consensus += aln['ref_seq'][ref_pos:ref_pos+l]
                ref_pos += l
            elif t=='I':
                for i in range(l):
                    R[i+consensus_pos] = '-'
                consensus += aln['query_seq'][query_pos:query_pos+l]
                query_pos += l

        new_block.sequences[aln['ref_name']] = R
        new_block.sequences[aln['query_name']] = S
        new_block.consensus = np.array(list(consensus))

        return new_block

    @classmethod
    def concatenate(cls, b1, b2):
        new_block = cls()
        assert b1.sequences.keys()==b2.sequences.keys()
        l1 = len(b1)
        new_block.consensus = np.concatenate((b1.consensus, b2.consensus))
        new_block.sequences = {s:dict(x) for s,x in b1.sequences.items()}
        for s in new_block.sequences:
            new_block.sequences[s].update({p+l1:c for p,c in b2.sequences[s].items()})

        return new_block


    @classmethod
    def from_cluster_alignment(cls, aln):
        new_block = cls()
        consensus = ""
        new_block.sequences = {}
        ref_pos = aln['ref_start']
        query_pos = 0
        consensus_pos = 0
        ref_map = [(ref_pos, 0)]
        query_map = [(query_pos, 0)]

        # additional changes
        R = {}
        S = {}
        tmpCigar = Cigar(aln["cigar"])

        for l,t in tmpCigar.items():
            if t in ['S', 'H']:
                query_pos+=l
            elif t=='M':
                ref_block = np.array(list(aln['ref_seq'][ref_pos:ref_pos+l]))
                query_block = np.array(list(aln['query_seq'][query_pos:query_pos+l]))
                diff = np.where(np.array(ref_block!=query_block))[0]
                for i in diff:
                    S[i+consensus_pos] = query_block[i]
                consensus+=aln['ref_seq'][ref_pos:ref_pos+l]
                ref_pos += l
                query_pos += l
            elif t=='D':
                for i in range(l):
                    S[i+consensus_pos] = '-'
                consensus+=aln['ref_seq'][ref_pos:ref_pos+l]
                ref_pos += l
            elif t=='I':
                for i in range(l):
                    R[i+consensus_pos] = '-'
                consensus+=aln['query_seq'][query_pos:query_pos+l]
                query_pos += l

            consensus_pos = len(consensus)

            ref_map.append((ref_pos, consensus_pos-ref_pos))
            query_map.append((query_pos, consensus_pos-query_pos))

        ref_map = np.array(ref_map).T
        query_map = np.array(query_map).T

        for transform, extra_mods, sub_cluster in [(ref_map, R, aln['ref_cluster']),
                                                   (query_map, S, aln['query_cluster'])]:
            for s,muts in sub_cluster.items():
                old_pos = np.array(list(muts.keys()))
                new_pos = old_pos + transform[1][np.searchsorted(transform[0], old_pos, side='right')]
                tmp_seq = {newp:muts[oldp] for newp,oldp in zip(new_pos, old_pos)}
                for p,c in extra_mods.items():
                    if p in tmp_seq:
                        # if agrees with new consensus, remove as modification
                        if tmp_seq[p]==consensus[p]:
                            tmp_seq.pop(p)
                        #otherwise, don't do anything
                    else:
                        tmp_seq[p]=c
                new_block.sequences[s] = tmp_seq

        new_block.consensus = np.array(list(consensus))

        return new_block


    def extract(self, seqname, strip_gaps=True):
        tmp = np.copy(self.consensus)
        for p,s in self.sequences[seqname].items():
            tmp[p]=s

        if strip_gaps:
            return "".join(tmp[tmp!='-'])
        else:
            return "".join(tmp)


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
        n = len(self.consensus)-1
        for s in self.sequences:
            new_block.sequences[s] = {n-p: reverse_complement.get(c,c)
                                      for p,c in self.sequences[s].items()}
        return new_block


    def __getitem__(self, val):
        if isinstance(val, slice):
            b = Block()
            start = val.start or 0
            stop = val.stop or len(self.consensus)
            b.consensus = self.consensus[start:stop]
            for s, mods in self.sequences.items():
                b.sequences[s] = {p-start:c for p,c in self.sequences[s].items()
                                  if p>=start and p<stop}
            return b
        else:
            raise ValueError("item access not supported")

    def __len__(self):
        return len(self.consensus)


if __name__ == '__main__':
    seqs={'R1':"___ABCdEF123GHiJ", 'S1':'xxxABCD6789EFGHIJ',
          'R2':'ABCdEF123GHiJ', 'S2':'ABcdEFaaa23GHiJ'}
    aln = {'ref_seq':seqs["R1"], 'query_seq':seqs["S1"], 'cigar':'3S4M4I2M3D4M',
           'ref_name':'R1', 'query_name':'S1', 'ref_start':3}

    b = Block.from_alignment(aln)

    aln = {'ref_seq':seqs["R2"], 'query_seq':seqs["S2"], 'cigar':'6M3I1D6M',
            'ref_name':'R2', 'query_name':'S2', 'ref_start':0}

    c = Block.from_alignment(aln)

    caln = {"ref_seq":"".join(b.consensus), "query_seq":"".join(c.consensus),
            "ref_cluster":b.sequences, "query_cluster":c.sequences,
            "ref_start":0, "cigar":"4M4D2M3I7M"}

    d  = Block.from_cluster_alignment(caln)

    for s in b.sequences:
        print(s, b.extract(s), seqs[s], b.extract(s)==seqs[s].strip('_x'))

    for s in c.sequences:
        print(s, c.extract(s), seqs[s], c.extract(s)==seqs[s].strip('_x'))

    for s in d.sequences:
        print(s, d.extract(s), seqs[s], d.extract(s)==seqs[s].strip('_x'))
