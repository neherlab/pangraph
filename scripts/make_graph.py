import os
from Bio import SeqIO, Phylo
from graph import Graph
from util import parse_paf
from betatree import betatree

#cluster = 'cluster_065_n_12' #'cluster_093_n_5'
#cluster = 'cluster_000_n_157' #'cluster_093_n_5'
cluster = 'cluster_022_n_19'
working_dir = cluster+'_dir'
if not os.path.isdir(working_dir):
	os.mkdir(working_dir)

seqs = list(SeqIO.parse(f'inc_clusters/{cluster}.fasta','fasta'))
unique_seqs = {seq.id:seq for seq in seqs}
seqs = list(unique_seqs.values())
bT = betatree.betatree(len(seqs))
bT.coalesce()
T = bT.BioTree

for (n,seq) in zip(T.get_terminals(), seqs):
	n.graph = Graph.from_sequence(seq.id, str(seq.seq).upper())
	n.name = seq.id #n.graph.sequences[seq.id][0][0]
	n.fasta_fname = os.path.join(*[working_dir, n.name+'.fasta'])
	n.graph.to_fasta(n.fasta_fname)

node_count = 0
for n in T.get_nonterminals(order='postorder'):
	n.name = f'NODE_{node_count:07d}'
	n.graph = Graph.fuse(n.clades[0].graph, n.clades[1].graph)
	os.system(f"minimap2 -x asm5 -D -c  {n.clades[0].fasta_fname} {n.clades[1].fasta_fname} > {working_dir}/{n.name}.paf")
	paf = parse_paf(f'{working_dir}/{n.name}.paf')
	paf.sort(key=lambda x:-x['aligned_bases'])
	merged_blocks = set()
	for hit in paf:
		if hit['query']['name'] in merged_blocks \
		   or hit['ref']['name'] in merged_blocks \
		   or hit['ref']['name']==hit['query']['name']:
			continue
		n.graph.merge_hit(hit)
		merged_blocks.add(hit['ref']['name'])
		merged_blocks.add(hit['query']['name'])

	n.fasta_fname = os.path.join(*[working_dir, n.name+'.fasta'])
	n.graph.to_fasta(n.fasta_fname)

Phylo.draw_ascii(T)
for (n,seq) in zip(T.get_terminals(), seqs):
	print(n.name)
	orig = str(seq.seq).upper()
	rec = T.root.graph.extract(n.name)
	print(orig==rec)
	if (orig!=rec):
		for i in range(len(orig)//100):
			if (orig[i*100:(i+1)*100]!=rec[i*100:(i+1)*100]):
				# for j,o,r in zip(range(100), orig[i*100:(i+1)*100], rec[i*100:(i+1)*100]):
				# 	if o!=r:
				# 		print(i*100+j)

				print(i,orig[i*100:(i+1)*100])
				print(i,rec[i*100:(i+1)*100])
		import ipdb; ipdb.set_trace()
	print(len(T.root.graph.extract(n.name)))

T.root.graph.to_json(os.path.join(working_dir, 'graph.json'))
