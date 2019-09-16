import argparse, gzip
from collections import defaultdict
from enum import IntEnum

import numpy as np
from matplotlib import pyplot as plt

# Rows of the PAF table
class PAF(IntEnum):
    qname  = 0
    qlen   = 1
    qbeg   = 2
    qend   = 3
    strand = 4
    rname  = 5
    rlen   = 6
    rbeg   = 7
    rend   = 8
    nmatch = 9
    alen   = 10
    mapq   = 11

class Hit():
    """docstring for Hit"""
    def __init__(self, col):
        self.qry    = (col[PAF.qname], int(col[PAF.qbeg]), int(col[PAF.qend]))
        self.ref    = (col[PAF.rname], int(col[PAF.rbeg]), int(col[PAF.rend]))
        self.nmatch = int(col[PAF.nmatch])
        self.alen   = int(col[PAF.alen])
        self.strand = col[PAF.strand]

    @property
    def qname(self):
        return self.qry[0]

    @property
    def rname(self):
        return self.ref[0]

    @property
    def aln(self):
        return (self.nmatch, self.alen, self.strand)

class Cluster(object):
    """docstring for Cluster"""
    def __init__(self, hit):
        super(Cluster, self).__init__()
        self.foot_print = [ 0, hit.alen ]
        self.fragments  = { hit.qname : { 'in_seq' : (hit.qry[1], hit.qry[2], hit.alen), 'in_cluster' : (0, hit.alen), 'strand' : 1 }, \
                            hit.rname : { 'in_seq' : (hit.ref[1], hit.ref[2], hit.alen), 'in_cluster' : (0, hit.alen), 'strand' : 1 } }

        if hit.strand == '-':
            self.fragments[hit.rname]['strand'] = -1

    def check_overlap(self, hit):
        if hit.qname in self.fragments or hit.rname in self.fragments:
            # sort partners by number of aligned bases
            # in most cases there will be just 1 -> no sorting necessary
            if hit.qname in self.fragments and hit.rname in self.fragments:
                refhit, qryhit = sorted([hit.ref, hit.qry], key = lambda x : self.fragments[x[0]]["in_seq"][2])
            elif hit.qryname in self.fragments:
                refhit, qryhit = hit.qry, hit.ref
            else:
                qryhit, refhit = hit.ref, hit.qry
            aln_info = hit.aln

            sh, bh, eh = refhit
            bc, ec, _a = self.fragments[sh]["in_seq"]
            if (eh < bc or bh > ec): # no overlap
                overlap = None
            elif (bh >= bc and eh <= ec):
                overlap = "contained"
            elif (bh < bc and eh > ec):
                overlap = "super"
            elif (bh < bc and eh <= ec):
                overlap = "left_overhang"
            elif (bh >= bc and eh > ec):
                overlap = "right_overhang"
            else:
                raise NotImplementedError("This case is not foreseen")
            return (overlap, refhit, qryhit, aln_info)

        else: # no overlap.
            return None

    def convert_coordinates(self, ref, x):
        tmp   = self.fragments[ref]
        slope = (tmp['in_cluster'][1] - tmp['in_cluster'][0])/(tmp['in_seq'][1] - tmp['in_seq'][0])
        return tmp['in_cluster'][0] + int((x-tmp['in_seq'][0])*slope)

    def merge_hit(self, overlap, refhit, qryhit, aln_info):
        if overlap == 'contained':
            c1     = self.convert_coordinates(refhit[0], refhit[1])
            c2     = self.convert_coordinates(refhit[0], refhit[2])
            s1, s2 = qryhit[1], qryhit[2]

            if qryhit[0] in self.fragments:
                c1 = min(c1, self.fragments[qryhit[0]]['in_cluster'][0])
                c2 = max(c2, self.fragments[qryhit[0]]['in_cluster'][1])

            self.fragments[qryhit[0]] = { 'in_seq' : (s1, s2, aln_info[2]), 'in_cluster' : (c1, c2) }
            print(refhit, self.fragments[refhit[0]])

def new_cluster(hit):
    c = {}

# ------------------------------------------------------------------------
# Main point of entry

def openany(fname, mode = 'r'):
    if fname.endswith('.gz'):
        return gzip.open(fname, mode)
    else:
        return open(fname, mode)

def main(args):
    maps = defaultdict(list)
    cov  = {}
    seq_names = {}

    i = 0
    with openany(args.fasta) as fh:
        for line in fh:
            if line.startswith('>'):
                seq_names[line[1:].strip().split()[0]] = i
                i += 1

    hits = []
    with openany(args.paf) as fh:
        for line in fh:
            entry   = line.strip().split()
            i1, i2  = seq_names[entry[PAF.qname]], seq_names[entry[PAF.rname]]
            if i1 < i2:
                hits.append(Hit(entry))
                maps[i1].append((int(entry[PAF.qbeg]), int(entry[PAF.qend])))
                maps[i2].append((int(entry[PAF.rbeg]), int(entry[PAF.rend])))

    for s in maps:
        maps[s].sort()
        tmp    = np.array(sorted([(x[0], 1) for x in maps[s]] + [(x[1], -1) for x in maps[s]]))
        cov[s] = (tmp[:,0], np.cumsum(tmp[:,1]))

    # sort by number of aligned bases
    hits.sort(key = lambda x: x.nmatch)

    tophit = hits[0]
    # a, b = tophit[0][0], tophit[0][1]

    focal_cluster = Cluster(tophit)
    print(focal_cluster.fragments)
    for p in hits[1:]:
        o = focal_cluster.check_overlap(p)
        if o:
            print(o)
            if o[2][0] == 4:
                import ipdb; ipdb.set_trace()
            if o[0] =='contained':
                focal_cluster.merge_hit(*o)

    # clusters = []
    # pair_to_clusters = {}
    # min_length = 100

    # for hit in hits:
    #   relevant_clusters = [c for c in clusters
    #                        if hit[0][0] in c.seqs or hit[1][0] in c.seqs]

    #   overlaps = []
    #   for c in relevant_clusters:
    #       overlaps.append((c,get_overlap(c, hit, min_length = min_length))

    #   if len(overlaps):
    #       overlaps.sort(key=lambda x:-x[-1])
    #       determine_merge(overlaps, c)
    #   else:
    #       create(new_cluster)

    #   # for c in relevant_clusters:


parser = argparse.ArgumentParser(description = "", usage = "parse and cluster a paf file")
parser.add_argument("fasta", type = str, help = "fasta file")
parser.add_argument("paf", type = str, help = "paf file")

if __name__ == '__main__':
    args = parser.parse_args()
    main(args)
