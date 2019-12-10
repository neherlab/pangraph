import numpy as np

from cigar import Cigar
from util  import partition_cigar

# ------------------------------------------------------------------------
# Global variables

reverse_complement = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}

def genblkid():
    alphabet = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    return "".join([alphabet[i] for i in np.random.randint(0, len(alphabet), 10)])

def mkarray(x):
    return np.array(list(x))

def mkblk(aln):
    R = {}
    Q = {}

    blkseq = ""
    refpos = 0 # aln['ref_start']
    qrypos = 0 # aln['qry_start']
    blkpos = 0

    refmap = [(refpos, 0)]
    qrymap = [(qrypos, 0)]

    refseq = aln['ref_seq']
    qryseq = aln['qry_seq']

    cigar = Cigar(aln["cigar"])

    for l, t in cigar.items():
        if t in ['S', 'H']:
            qrypos += l

        elif t == 'M':
            refblk = mkarray(refseq[refpos:refpos+l])
            qryblk = mkarray(qryseq[qrypos:qrypos+l])
            diff   = np.where(np.array(refblk != qryblk))[0]
            for i in diff:
                Q[i+blkpos] = qryblk[i]

            blkseq += aln['ref_seq'][refpos:refpos+l]
            refpos += l
            qrypos += l

        elif t == 'D':
            for i in range(l):
                Q[i+blkpos] = '-'
            blkseq += aln['ref_seq'][refpos:refpos+l]
            refpos += l

        elif t == 'I':
            for i in range(l):
                R[i+blkpos] = '-'
            blkseq += aln['qry_seq'][qrypos:qrypos+l]
            qrypos += l

        blkpos = len(blkseq)

        refmap.append((refpos, blkpos-refpos))
        qrymap.append((qrypos, blkpos-qrypos))

    refmap = np.array(refmap).T
    qrymap = np.array(qrymap).T

    return R, Q, blkseq, refmap, qrymap


# def mkblk(aln, cutoff=500):
#     R = {}
#     Q = {}

#     blkseq = ""
#     ql, qr = 0, 0
#     rl, rr = 0, 0

#     refpos = 0 # aln['ref_start']
#     qrypos = 0 # aln['qry_start']
#     blkpos = 0

#     refmap = [(refpos, 0)]
#     qrymap = [(qrypos, 0)]

#     refseq = aln['ref_seq']
#     qryseq = aln['qry_seq']

#     refblks = []
#     qryblks = []

#     print("######", partition_cigar(aln["cigar"]))
#     cigar = Cigar(aln["cigar"])

#     for l, t in cigar.items():
#         if t in ['S', 'H']:
#             if l >= cutoff:
#                 ql = qr
#             # qrypos += l
#             qr += l

#         elif t == 'M':
#             refblk = mkarray(refseq[rr:rr+l])
#             qryblk = mkarray(qryseq[qr:qr+l])
#             diff   = np.where(np.array(refblk != qryblk))[0]
#             for i in diff:
#                 Q[i+qr] = qryblk[i]

#             blkseq += aln['ref_seq'][rr:rr+l]
#             rr += l
#             qr += l

#         elif t == 'D':
#             if l >= cutoff:
#                 # Need to make a new block
#             else:
#                 for i in range(l):
#                     Q[i+blkpos] = '-'
#                 blkseq += aln['ref_seq'][rr:rr+l]
#             rr += l

#         elif t == 'I':
#             if l >= cutoff:
#             else:
#                 for i in range(l):
#                     R[i+blkpos] = '-'
#                 blkseq += aln['qry_seq'][qr:qr+l]
#                 qr += l

#         br = len(blkseq)

#         refmap.append((rr, br-rr))
#         qrymap.append((qr, br-qr))

#     refmap = np.array(refmap).T
#     qrymap = np.array(qrymap).T

#     return R, Q, blkseq, refmap, qrymap

# ------------------------------------------------------------------------
# Block class

class Block(object):
    """docstring for Block"""
    def __init__(self):
        super(Block, self).__init__()
        self.muts = {}
        self.seq  = None
        self.id   = genblkid()

    @classmethod
    def from_seq(cls, name, seq):
        new_blk      = cls()
        new_blk.seq  = np.array(list(seq))
        new_blk.muts = {name:{}}
        return new_blk

    @classmethod
    def cat(cls, blks):
        new_blk = cls()
        assert all([blks[0].muts.keys()==b2.muts.keys() for b2 in blks[1:]])

        new_blk.seq = np.concatenate([b.seq for b in blks])
        new_blk.muts = {s:dict(x) for s,x in blks[0].muts.items()}
        offset = len(blks[0])
        for b in blks[1:]:
            for s in new_blk.muts:
                new_blk.muts[s].update({p+offset:c for p,c in b.muts[s].items()})
            offset += len(b)

        return new_blk

    @classmethod
    def from_aln(cls, aln, verify=True):
        # print("###### PARTITION CIGAR ########")
        # print("CIGAR:", aln['cigar'])
        qrys, refs, blks = partition_cigar(aln['cigar'], aln['qry_seq'], aln['ref_seq'])
        # print("Intervals", qrys, refs)

        if verify:
            # assert len(refs) == 1 and len(qrys) == 1 and len(blks) == 1, "non unity number of blocks"
            cseq = "".join(n for b in blks for n in b[0])
            R, Q, blkseq, refmap, qrymap = mkblk(aln)
            # print(f"OUR SEQ: {cseq[0:100]}")
            # print(f"RCH SEQ: {blkseq[0:100]}")
            assert blkseq == cseq, "unequal seqs"
            if len(blks) == 1:
                assert np.all(qrymap == blks[0][1][1]), "unequal qrymap"
                assert np.all(refmap == blks[0][2][1]), "unequal refmap"
            else:
                pass
                # for b in blks:
                    # print("--> NQ:", b[1][0])
                    # print("--> OQ:", Q)
                    # print("--> NR:", b[2][0])
                    # print("--> OR:", R)
                    # import ipdb
                    # ipdb.set_trace()

            # cseq = "".join(n for n in blks[0][0])
            # print("TYPE | OURS | RICHARD")
            # print(f"--> REFMAP: {refmap} | {blks[0][2][1]}")
            # print(f"--> QRYMAP: {qrymap} | {blks[0][1][1]}")
            # print(f"OUR SEQ: {cseq[0:100]}")
            # print(f"RCH SEQ: {blkseq[0:100]}")
            # assert blkseq == cseq, "unequal seqs"

        def updatemuts(blk, xtramuts, xmap, omuts, ival):
            seq = blk.seq
            # Iterate over all sequences in the block
            for iso, muts in omuts.items():

                # Shift position by indel #'s
                opos = mkarray(m for m in muts.keys() if m >= ival[0] and m < ival[1])
                # print("----> XMAP:", xmap)
                # print("----> OPOS:", opos)
                # print("----> SEARCH:", np.searchsorted(xmap[0], opos, side='right'))
                npos = opos + xmap[1][np.searchsorted(xmap[0], opos, side='right')]

                assert (np.max(npos) if len(npos) > 0 else 0) < len(seq)
                newmuts = {np : muts[op] for np, op in zip(npos, opos)}
                # print("----> NEWMUTS:", xtramuts)
                for p, n in xtramuts.items():
                    if p in newmuts:
                        if newmuts[p] == seq[p]:
                            newmuts.pop(p)
                    else:
                        newmuts[p] = n
                blk.muts[iso] = newmuts
            return blk

        newblks = []
        for i, blk in enumerate(blks):
            newblk      = cls()
            newblk.muts = {}

            newblk.seq, qry, ref = blk
            Q, qrymap = qry
            R, refmap = ref
            # print("====== MUTATION TRANSFER ======")
            # print("--> Intervals", qrys[i], ":", refs[i])
            # print("--> Q", Q)
            # print("--> R", R)
            # print("--> QRYMUT", aln["qry_cluster"])
            # print("--> REFMUT", aln["ref_cluster"])

            if qrys[i] is not None:
                # print("--> ADDING Q")
                newblk = updatemuts(newblk, Q, qrymap, aln['qry_cluster'], qrys[i])
            if refs[i] is not None:
                # print("--> ADDING R")
                newblk = updatemuts(newblk, R, refmap, aln['ref_cluster'], refs[i])
            newblks.append(newblk)

            # print("NEWMUT", newblk.muts)

        qryblks = [nb for i, nb in enumerate(newblks) if qrys[i] is not None]
        if aln['orientation'] == -1:
            qryblks = qryblks[::-1]
        refblks = [nb for i, nb in enumerate(newblks) if refs[i] is not None]
        return newblks, qryblks, refblks

    def extract(self, iso, strip_gaps=True, verbose=False):
        tmp = np.copy(self.seq)
        # if verbose:
        #     print("----> ISO", iso)
        #     print("----> MUTS", self.muts.keys())
        for p, s in self.muts[iso].items():
            tmp[p] = s

        assert len(tmp) > 0, "empty sequence"
        if strip_gaps:
            return "".join(tmp[tmp != '-'])
        else:
            return "".join(tmp)

    def copy(self, keep_name=False):
        b = Block()
        if keep_name:
            b.name = self.name
        b.seq = np.copy(self.seq)
        b.muts = {s:dict(self.muts[s]) for s in self.muts}
        return b

    def rev_cmpl(self):
        from Bio import Seq

        new_blk = Block()
        new_blk.seq = np.array(list(Seq.reverse_complement("".join(self.seq))))
        n = len(self.seq)-1
        for s in self.muts:
            new_blk.muts[s] = {n-p: reverse_complement.get(c,c) for p,c in self.muts[s].items()}
        return new_blk


    def __getitem__(self, val):
        if isinstance(val, slice):
            b = Block()
            start = val.start or 0
            stop = val.stop or len(self.seq)
            b.seq = self.seq[start:stop]
            for s, _ in self.muts.items():
                b.muts[s] = {p-start:c for p,c in self.muts[s].items() if p>=start and p<stop}
            return b
        else:
            raise ValueError("item access not supported")

    def __len__(self):
        return len(self.seq)

if __name__ == '__main__':
    # seqs={'R1':"___ABCdEF123GHiJ",\
    #       'S1':"___ABCDEFGHi"}
    # aln = {'ref_seq':seqs["R1"], 'qry_seq':seqs["S1"], 'cigar':'9M3D3M1D',
    #        'ref_name':'R1', 'qry_name':'S1', 'ref_start':0}

    # b = Block.from_aln(aln)

    # for s in b.sequences:
        # print(s, b.extract(s), seqs[s], b.extract(s)==seqs[s].strip('_x'))

    # seqs={'R1':"___ABCdEF123GHiJ",\
    #       'S1':"___ABCDEFGHiKLM"}
    # aln = {'ref_seq':seqs["R1"], 'qry_seq':seqs["S1"], 'cigar':'9M3D3M1D3I',
    #        'ref_name':'R1', 'qry_name':'S1', 'ref_start':0}

    # c = Block.from_aln(aln)

    # for s in c.sequences:
    #     print(s, b.extract(s), seqs[s], b.extract(s)==seqs[s].strip('_x'))

    seqs={'R1':"___ABCdEF123GHiJ",\
          'S1':"___ABCDefGHiKLM"}
    aln = {'ref_seq':seqs["R1"], 'qry_seq':seqs["S1"], 'cigar':'9M3D3M1D3I',
           'ref_name':'R1', 'qry_name':'S1', 'ref_start':0}

    c = Block.from_aln(aln)

    for s in c.sequences:
        print(s, c.extract(s), seqs[s], c.extract(s)==seqs[s].strip('_x'))
