import numpy as np
import numpy.random as rng

from utils import parsecigar, wcpair, asarray

# ------------------------------------------------------------------------
# Helper functions

rng.seed(0) # For deterministic block names
def randomid():
    alphabet = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    return "".join([alphabet[i] for i in rng.randint(0, len(alphabet), 10)])

# ------------------------------------------------------------------------
# Block class

class Block(object):
    """docstring for Block"""

    def __init__(self):
        super(Block, self).__init__()
        self.id   = randomid()
        self.seq  = None
        self.muts = {}

    # --- Class methods ---

    @classmethod
    def fromseq(cls, name, seq):
        new_blk      = cls()
        new_blk.seq  = asarray(seq)
        new_blk.muts = {name:{}}

        return new_blk

    @classmethod
    def cat(cls, blks):
        nblk = cls()
        assert all([blks[0].muts.keys()==b2.muts.keys() for b2 in blks[1:]])

        nblk.seq  = np.concatenate([b.seq for b in blks])
        nblk.muts = { s:dict(x) for s,x in blks[0].muts.items() }
        offset    = len(blks[0])
        for b in blks[1:]:
            for s in nblk.muts:
                nblk.muts[s].update({p+offset:c for p,c in b.muts[s].items()})
            offset += len(b)

        return nblk

    @classmethod
    def fromaln(cls, aln):
        def updatemuts(blk, xtramuts, xmap, omuts, ival):
            seq = blk.seq
            # Iterate over all sequences in the block
            for iso, muts in omuts.items():

                # Shift position by indel #'s
                opos = asarray(m for m in muts.keys() if m >= ival[0] and m < ival[1])
                npos = opos + xmap[1][np.searchsorted(xmap[0], opos, side='right')]

                assert (np.max(npos) if len(npos) > 0 else 0) < len(seq)

                newmuts = {np : muts[op] for np, op in zip(npos, opos)}
                for p, n in xtramuts.items():
                    if p in newmuts:
                        if newmuts[p] == seq[p]:
                            newmuts.pop(p)
                    else:
                        newmuts[p] = n
                blk.muts[iso] = newmuts
            return blk

        qrys, refs, blks = parsecigar(aln['cigar'], aln['qry_seq'], aln['ref_seq'])

        newblks = []
        for i, blk in enumerate(blks):
            newblk      = cls()
            newblk.muts = {}

            newblk.seq, qry, ref = blk
            Q, qrymap = qry
            R, refmap = ref

            if qrys[i] is not None:
                newblk = updatemuts(newblk, Q, qrymap, aln['qry_cluster'], qrys[i])
            if refs[i] is not None:
                newblk = updatemuts(newblk, R, refmap, aln['ref_cluster'], refs[i])
            newblks.append(newblk)

        qryblks = [nb for i, nb in enumerate(newblks) if qrys[i] is not None]
        if aln['orientation'] == -1:
            qryblks = qryblks[::-1]
        refblks = [nb for i, nb in enumerate(newblks) if refs[i] is not None]

        return newblks, qryblks, refblks

    # --- Instance methods ---

    def copy(self, keepname=False):
        b = Block()
        if keepname:
            b.name = self.name
        b.seq  = np.copy(self.seq)
        b.muts = {s:dict(self.muts[s]) for s in self.muts}

        return b

    def extract(self, iso, strip_gaps=True, verbose=False):
        tmp = np.copy(self.seq)
        for p, s in self.muts[iso].items():
            tmp[p] = s

        assert len(tmp) > 0, "empty sequence"
        if strip_gaps:
            return "".join(tmp[tmp != '-'])
        else:
            return "".join(tmp)

    def revcmpl(self):
        from Bio import Seq

        nblk     = Block()
        nblk.seq = asarray(Seq.reverse_complement("".join(self.seq)))
        L        = len(self.seq)-1
        for s in self.muts:
            nblk.muts[s] = {L-p: wcpair.get(c,c) for p,c in self.muts[s].items()}

        return nblk

    # --- Python Operator Overloads ---

    def __len__(self):
        return len(self.seq)

    def __getitem__(self, val):
        if isinstance(val, slice):
            b     = Block()
            start = val.start or 0
            stop  = val.stop or len(self.seq)
            b.seq = self.seq[start:stop]
            for s, _ in self.muts.items():
                b.muts[s] = {p-start:c for p,c in self.muts[s].items() if p>=start and p<stop}
            return b
        else:
            raise ValueError("item access not supported")
