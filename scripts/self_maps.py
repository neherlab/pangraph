import os

import numpy as np

import Cigar
import graph

# -----------------------------------------------------------------------------
# Functions

def self_maps(G, qry, ref, out):
    os.system(f"minimap2 -x asm5 -D -c {qry} {ref} 1> {out}.paf 2>log")

    paf = parse_paf(f"{out}.paf")
    paf.sort(key=lambda x:-x['aligned_bases'])

    merged_blocks = set()
    for hit in paf:
        graph.merge_hit(hit)
        merged_blocks.add(hit['ref']['name'])
        merged_blocks.add(hit['qry']['name'])

    G.prune_empty()
    return G

# -----------------------------------------------------------------------------
# Main point of entry

if __name__ == "__main__":
