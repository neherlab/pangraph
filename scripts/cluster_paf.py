import argparse, gzip
import numpy as np
from collections import defaultdict
from matplotlib import pyplot as plt

def myopen(fname, mode='r'):
    if fname.endswith('.gz'):
        return gzip.open(fname, mode)
    else:
        return open(fname, mode)


class Cluster(object):
    """docstring for Cluster"""
    def __init__(self, hit):
        super(Cluster, self).__init__()
        self.foot_print = [0,np.max([h[2]-h[1] for h in hit[:2]])]
        self.fragments = {h[0]:{'in_seq':(h[1],h[2], hit[2][1]), 'in_cluster':(0,self.foot_print[-1]), 'strand':1}
                          for h in hit[:2]}
        if hit[2][-1]=='-':
            self.fragments[hit[1][0]]['strand']=-1


    def check_overlap(self,hit):
        if hit[0][0] in self.fragments or hit[1][0] in self.fragments:
            # sort partners by number of aligned bases
            # in most cases there will be just 1 -> no sorting necessary
            if hit[0][0] in self.fragments and hit[1][0] in self.fragments:
                reference_hit, other = sorted([h for h in hit[:2] if h[0] in self.fragments],
                                    key=lambda x:self.fragments[x[0]]["in_seq"][2])
            elif hit[0][0] in self.fragments:
                reference_hit, other = hit[:2]
            else:
                other, reference_hit = hit[:2]
            aln_info = hit[2]

            sh,bh,eh = reference_hit
            bc,ec,_a = self.fragments[sh]["in_seq"]
            if (eh<bc or bh>ec): # no overlap
                overlap = None
            elif (bh>=bc and eh<=ec):
                overlap = "contained"
            elif (bh<bc and eh>ec):
                overlap = "super"
            elif (bh<bc and eh<=ec):
                overlap = "left_overhang"
            elif (bh>=bc and eh>ec):
                overlap = "right_overhang"
            else:
                raise NotImplementedError("This case is not foreseen")
            return (overlap, reference_hit, other, aln_info)

        else: # no overlap.
            return None

    def convert_coordinates(self,ref, x):
        tmp = self.fragments[ref]
        slope = (tmp['in_cluster'][1]-tmp['in_cluster'][0])/(tmp['in_seq'][1]-tmp['in_seq'][0])
        # print(slope, tmp, x, tmp['in_cluster'][0] + int((x-tmp['in_seq'][0])*slope))
        return tmp['in_cluster'][0] + int((x-tmp['in_seq'][0])*slope)

    def merge_hit(self, overlap, reference_hit, other, aln_info):
        if overlap=='contained':
            c1 = self.convert_coordinates(reference_hit[0], reference_hit[1])
            c2 = self.convert_coordinates(reference_hit[0], reference_hit[2])
            s1 = other[1]
            s2 = other[2]
            if other[0] in self.fragments:
                c1 = min(c1, self.fragments[other[0]]['in_cluster'][0])
                c2 = max(c2, self.fragments[other[0]]['in_cluster'][1])
            self.fragments[other[0]] = {'in_seq': (s1,s2, aln_info[2]),
                                        'in_cluster':(c1,c2)}
            print(reference_hit, self.fragments[reference_hit[0]])

def new_cluster(hit):
    c = {}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = "", usage="parse and cluster a paf file")
    parser.add_argument("fasta", type=str, help="fasta file")
    parser.add_argument("paf", type=str, help="paf file")
    args=parser.parse_args()

    mappings = defaultdict(list)
    coverage = {}
    seq_names = {}
    i=0
    with open(args.fasta) as fh:
        for line in fh:
            if line[0]=='>':
                seq_names[line[1:].strip().split()[0]] = i
                i+=1

    hits = []
    with myopen(args.paf) as fh:
        for line in fh:
            entries = line.strip().split()
            i1,i2 = seq_names[entries[0]], seq_names[entries[5]]
            if i1<i2:
                hits.append([(i1, int(entries[2]), int(entries[3])),
                             (i2, int(entries[7]), int(entries[8])),
                             (int(entries[9]), int(entries[10]), entries[4]=='+')])
                mappings[i1].append((int(entries[2]), int(entries[3])))
                mappings[i2].append((int(entries[7]), int(entries[8])))


    for s in mappings:
        mappings[s].sort()
        tmp = np.array(sorted([(x[0],1) for x in mappings[s]] + [(x[1],-1) for x in mappings[s]]))
        coverage[s] = (tmp[:,0], np.cumsum(tmp[:,1]))
        # plt.plot(coverage[s][0], coverage[s][1])

    # sort by number of aligned bases
    hits.sort(key = lambda x:-x[2][0])

    focus_hit = hits[0]
    a, b = focus_hit[0][0], focus_hit[0][1]

    focal_cluster = Cluster(focus_hit)
    print(focal_cluster.fragments)
    for p in hits[1:]:
        o = focal_cluster.check_overlap(p)
        if o:
            print(o)
            if o[2][0]==4:
                import ipdb; ipdb.set_trace()
            if o[0]=='contained':
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
