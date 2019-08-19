import argparse, gzip
import numpy as np
from collections import defaultdict

def myopen(fname, mode='r'):
    if fname.endswith('.gz'):
        return gzip.open(fname, mode)
    else:
        return open(fname, mode)


def parse_paf(fname):
    hits = []
    with myopen(fname) as fh:
        for line in fh:
            entries = line.strip().split()
            hit = {'query':{'name':entries[0], 'start':int(entries[2]), 'end':int(entries[3])},
                   'ref':  {'name':entries[5], 'start':int(entries[7]), 'end':int(entries[8])},
                   'aligned_bases':int(entries[9]), 'aligned_length':int(entries[10]),
                   'orientation':1 if entries[4]=='+' else -1, 'mapping_quality':int(entries[11])}
            for extra in entries[12:]:
                if extra.startswith('cg:'):
                    hit['cigar'] = extra.split(':')[-1]
            hits.append(hit)
    return hits
