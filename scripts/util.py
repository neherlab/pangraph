import argparse, gzip
import csv
import numpy as np
from collections import defaultdict

def myopen(fname, mode='r'):
    if fname.endswith('.gz'):
        return gzip.open(fname, mode)
    else:
        return open(fname, mode)

def parse_paf(fname):
    assert fname.endswith(".paf")
    hits = []
    with myopen(fname) as fh:
        for line in fh:
            row = line.strip().split()
            hit = {'qry': {'name' : row[0],
                           'start' : int(row[2]),
                           'end' : int(row[3])},
                   'ref': {'name' : row[5],
                           'start' : int(row[7]),
                           'end' : int(row[8])},
                   'aligned_bases' : int(row[9]),
                   'aligned_length' : int(row[10]),
                   'mapping_quality':int(row[11]),
                   'orientation':1 if row[4]=='+' else -1}
            for extra in row[12:]:
                if extra.startswith('cg:'):
                    hit['cigar'] = extra.split(':')[-1]
            hits.append(hit)
    return hits

def parse_tsv(fname):
    assert fname.endswith(".tsv")
    hits = dict()
    with open(fname) as fh:
        rdr = csv.reader(fh, delimiter = "\t", quotechar = '"')
        for row in rdr:
            assert row[0] not in hits
            hits[row[0][:-6]] = (row[1], row[2])

    return dict(hits)

def parse_m8(fname, onlytophit=True):
    assert fname.endswith(".m8")
    hits = defaultdict(list)
    with open(fname) as fh:
        rdr = csv.reader(fh, delimiter = "\t", quotechar = '"')
        for row in rdr:
            hits[row[0]].append({ "hit" : row[1], "prcid" : float(row[2]), "eval" : float(row[11]), "beg" : int(row[6]), "end" : int(row[7]) })

        hits = dict(hits)
        if onlytophit:
            for qry in hits:
                hits[qry] = sorted(hits[qry], key = lambda x: x["eval"])[::-1]

    return hits

def parse_kmer(mtx):
    with open(mtx) as fh:
        nrows = int(fh.readline().strip())
        M = np.zeros((nrows, nrows), dtype=float)
        seq_names = []
        for li, line in enumerate(fh):
                e = line.strip().split()
                seq_names.append(e[0].split('/')[-1][:-3])
                M[li,:(li+1)] = [float(x) for x in e[1:]]

        Msym = M + M.T
        Msym -= np.eye(nrows)

    return Msym, np.array(seq_names)


def partition_cigar(aln, qryseq, refseq, cutoff=500):
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
