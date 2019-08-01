import os
from Bio import SeqIO

with open('../../GROUP/data/plasmids/kmer_groups/inc_group_clusters.csv') as fh:
    for li, line in enumerate(fh):
        files = line.strip().split(',')
        out_file = f'inc_clusters/cluster_{li:03d}_n_{len(files)}.fasta'
        with open(out_file, 'w') as outfh:
            for fname in files:
                SeqIO.write(SeqIO.parse(fname, 'fasta'), outfh, 'fasta')

