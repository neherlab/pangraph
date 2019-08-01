from Bio import SeqIO
from graph import Graph
from util import parse_paf


seqs = list(SeqIO.parse('inc_clusters/cluster_067_n_2.fasta','fasta'))

graphs = [Graph.from_sequence(s.id, str(s.seq).upper()) for s in seqs]

g = Graph.fuse(graphs[0], graphs[1])
paf = parse_paf('mapping/cluster_067_n_2.paf')

paf[0]['ref']['name'] = g.sequences[paf[0]['ref']['name']][0][0]
paf[0]['query']['name'] = g.sequences[paf[0]['query']['name']][0][0]
g.merge_hit(paf[0])


