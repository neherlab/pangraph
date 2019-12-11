import gzip
import numpy as np

from enum import IntEnum

# ------------------------------------------------------------------------
# Global Enums

class Strand(IntEnum):
    Plus  = +1
    Minus = -1

wcpair = {'A' : 'T', 'T': 'A', 'C' : 'G', 'G': 'C'}

# ------------------------------------------------------------------------
# Global Functions

def panic(msg):
    raise SystemExit(f"Panic: {msg}")

def asarray(x):
    return np.array(list(x))

# ------------------------------------------------------------------------
# Helper Functions

def openany(path, mode='r'):
    if path.endswith('.gz'):
        return gzip.open(path, mode)
    else:
        return open(path, mode)

def parsepaf(path):
    assert path.endswith(".paf")
    hits = []
    with openany(path) as fh:
        for line in fh:
            row = line.strip().split()
            hit = {'qry': {'name'    : row[0],
                           'start'   : int(row[2]),
                           'end'     : int(row[3])},
                   'ref': {'name'    : row[5],
                           'start'   : int(row[7]),
                           'end'     : int(row[8])},
                   'aligned_bases'   : int(row[9]),
                   'aligned_length'  : int(row[10]),
                   'mapping_quality' : int(row[11]),
                   'orientation'     : 1 if row[4]=='+' else -1}
            for xtra in row[12:]:
                if xtra.startswith('cg:'):
                    hit['cigar'] = xtra.split(':')[-1]
            hits.append(hit)
    return hits

def parsecigar(aln, qryseq, refseq, cutoff=500):
    from cigar import Cigar

    aln = Cigar(aln)

    lq, rq = 0, 0
    lr, rr = 0, 0
    refs = []
    qrys = []
    blks = []

    R, Q = {}, {}
    blkseq = ""
    blkpos = 0
    refmap = [(rr, blkpos-rr)]
    qrymap = [(rq, blkpos-rq)]

    def push(qval=None, rval=None):
        nonlocal R, Q, blkseq, blkpos, refmap, qrymap
        assert not (qval is None and rval is None)

        def f(xs, x):
            if x is None:
                xs.append(None)
                return True
            else:
                l, r = zip(x)
                if l < r:
                    xs.append(x)
                    return True
                return False

        hasq = f(qrys, qval)
        hasr = f(refs, rval)

        if hasq or hasr:
            assert len(qrys) == len(refs)
            assert len(blkseq) > 0, "empty seq"
            blks.append((np.array(list(blkseq)), (Q, np.array(qrymap).T), (R, np.array(refmap).T)))

        R, Q = {}, {}
        blkseq = ""
        blkpos = 0
        refmap = [(rr, blkpos-rr)]
        qrymap = [(rq, blkpos-rq)]

    def recordbp():
        nonlocal blkpos, refmap, qrymap

        blkpos = len(blkseq)
        refmap.append((rr, blkpos-rr))
        qrymap.append((rq, blkpos-rq))

    for l, t in aln.items():
        if t in ['S', 'H']:
            if l >= cutoff:
                push((lq, rq), (lr, rr))

                blkseq = qryseq[rq:rq+l]
                rq += l
                recordbp()

                push((rq-l, rq), None)
                lq = rq
                lr = rr
            else:
                rq += l
                recordbp()

        elif t == 'M':
            rs = np.array(list(refseq[rr:rr+l]))
            qs = np.array(list(qryseq[rq:rq+l]))
            diff = np.where(np.array(rs != qs))[0]
            for i in diff:
                Q[i+blkpos] = qs[i]
            blkseq += refseq[rr:rr+l]

            rq += l
            rr += l

            recordbp()

        elif t == 'D':
            if l >= cutoff:
                push((lq, rq), (lr, rr))

                blkseq = refseq[rr:rr+l]
                rr += l
                recordbp()

                push(None, (rr-l, rr))
                lr = rr
                lq = rq
            else:
                for i in range(l):
                    Q[i+blkpos] = '-'
                blkseq += refseq[rr:rr+l]

                rr += l
                recordbp()

        elif t == 'I':
            if l >= cutoff:
                push((lq, rq), (lr, rr))

                blkseq = qryseq[rq:rq+l]
                rq += l
                recordbp()

                push((rq-l, rq), None)
                lq = rq
                lr = rr
            else:
                for i in range(l):
                    R[i+blkpos] = '-'
                blkseq += qryseq[rq:rq+l]

                rq += l
                recordbp()

    push((lq, rq), (lr, rr))

    assert len(qrys) == len(refs) and len(qrys) == len(blks)
    return qrys, refs, blks
