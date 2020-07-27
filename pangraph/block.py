import numpy as np
import numpy.random as rng

from collections import defaultdict
from .utils import parse_cigar, wcpair, as_array, as_string

# ------------------------------------------------------------------------
# Helper functions

RS = rng.RandomState(0)
def randomid():
    alphabet = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    randomid.N += 1
    name = "".join([alphabet[i] for i in RS.randint(0, len(alphabet), 10)])

    return name
randomid.N = 0

# ------------------------------------------------------------------------
# Block class

class Block(object):
    """docstring for Block"""

    def __init__(self, gen=True):
        super(Block, self).__init__()
        self.id   = randomid() if gen else 0
        self.seq  = None
        self.muts = {}

    # ------------------
    # properties

    @property
    def length(self):
        return len(self.seq)

    @property
    def depth(self):
        return len(self.muts)

    # ------------------
    # static methods

    @classmethod
    def from_seq(cls, name, seq):
        new_blk      = cls()
        new_blk.seq  = as_array(seq)
        new_blk.muts = {(name, 0):{}}

        return new_blk

    @classmethod
    def from_dict(cls, d):
        def unpack(key):
            t = tuple(key.split("?###?"))
            return (t[0], int(t[1]))

        B      = Block()
        B.id   = d['id']
        B.seq  = as_array(d['seq'])
        B.muts = {unpack(k):v for k, v in d['muts'].items()}

        return B

    @classmethod
    def cat(cls, blks):
        nblk = cls()
        assert all([blks[0].muts.keys() == b2.muts.keys() for b2 in blks[1:]])

        nblk.seq  = np.concatenate([b.seq for b in blks])
        nblk.muts = { s:dict(x) for s,x in blks[0].muts.items() }
        offset    = len(blks[0])
        for b in blks[1:]:
            for s in nblk.muts:
                nblk.muts[s].update({p+offset:c for p,c in b.muts[s].items()})
            offset += len(b)

        return nblk

    @classmethod
    def from_aln(cls, aln, debug=False):
        def updatemuts(blk, xtramuts, xmap, omuts, ival):
            seq = blk.seq
            # Iterate over all sequences in the block
            isomap = {}
            for iso, muts in omuts.items():
                # Shift position by indel #'s
                opos = as_array(m for m in muts.keys() if m >= ival[0] and m < ival[1])
                npos = opos + xmap[1][np.searchsorted(xmap[0], opos, side='right')]

                assert (np.max(npos) if len(npos) > 0 else 0) < len(seq)

                newmuts = {np : muts[op] for np, op in zip(npos, opos)}
                for p, n in xtramuts.items():
                    if p in newmuts:
                        if newmuts[p] == seq[p]:
                            newmuts.pop(p)
                    else:
                        newmuts[p] = n

                isomap[iso] = blk.push(iso, newmuts)

            return blk, isomap

        qrys, refs, blks = parse_cigar(aln['cigar'], aln['qry_seq'], aln['ref_seq'])

        # Iterate over all merged blocks and merge their sequences + mutations.
        newblks = []
        isomap  = defaultdict(lambda: {})
        for i, blk in enumerate(blks):
            newblk      = cls()
            newblk.muts = {}

            newblk.seq, qry, ref = blk
            Q, qrymap = qry
            R, refmap = ref

            if qrys[i] is not None:
                newblk, isomap[newblk.id][aln['qry_name']] = updatemuts(newblk, Q, qrymap, aln['qry_cluster'], qrys[i])
            if refs[i] is not None:
                newblk, isomap[newblk.id][aln['ref_name']] = updatemuts(newblk, R, refmap, aln['ref_cluster'], refs[i])

            newblks.append(newblk)

        isomap = dict(isomap)

        qryblks = [nb for i, nb in enumerate(newblks) if qrys[i] is not None]
        if aln['orientation'] == -1:
            qryblks = qryblks[::-1]
        refblks = [nb for i, nb in enumerate(newblks) if refs[i] is not None]

        return newblks, qryblks, refblks, isomap

    # --------------
    # methods

    def copy(self, keepname=False):
        b = Block()
        if keepname:
            b.name = self.name
        b.seq  = np.copy(self.seq)
        b.muts = {s:dict(self.muts[s]) for s in self.muts}

        return b

    def extract(self, iso, num, strip_gaps=True, verbose=False):
        tag = (iso, num)
        tmp = np.copy(self.seq)
        for p, s in self.muts[tag].items():
            tmp[p] = s

        assert len(tmp) > 0, "empty sequence"
        if strip_gaps:
            return as_string(tmp[tmp != '-'])
        else:
            return as_string(tmp)

    def len_of(self, iso, num):
        tag    = (iso, num)
        length = len(self.seq)
        gaplen = sum(1 for s in self.muts[tag].values() if s == '-')
        return length - gaplen

    def is_empty(self, iso, num, strip_gaps=True):
        tag = (iso, num)
        seq = np.copy(self.seq)
        # NOTE: This is a hack. Need to investigate the error that arises.
        if tag not in self.muts:
            return True

        for p, s in self.muts[tag].items():
            seq[p] = s

        return len(seq) == 0 or all(nuc == "-" for nuc in seq)

    def rev_cmpl(self):
        from Bio import Seq

        nblk     = Block()
        nblk.seq = as_array(Seq.reverse_complement("".join(self.seq)))
        L        = len(self.seq)-1
        for s in self.muts:
            nblk.muts[s] = {L-p: wcpair.get(c,c) for p,c in self.muts[s].items()}

        return nblk

    def push(self, iso, muts):
        tag = iso if isinstance(iso, tuple) else (iso, 0)

        while tag in self.muts:
            tag = (tag[0], tag[1]+1)

        self.muts[tag] = muts
        return tag

    def has(self, iso):
        # TODO: Removed assumption that isolates will be strictly ordered.
        # i.e. isolate num 0 can be removed but num 1 can still exist.
        # Test to see if this works
        for tag in self.muts.keys():
            if tag[0] == iso:
                return True
        return False

    def to_dict(self):
        def pack(key):
            return f"{key[0]}?###?{key[1]}"

        def fix(d):
            return {int(k):v for k, v in d.items()}

        return {'id'   : self.id,
                'seq'  : "".join(str(n) for n in self.seq),
                'muts' : {pack(k) : fix(v) for k, v in self.muts.items()}}

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
