import os
from Bio import SeqIO, Phylo
from graph import Graph
from util import parse_paf

cluster = 'cluster_065_n_12' #'cluster_093_n_5'
cluster = 'cluster_000_n_157' #'cluster_093_n_5'
#cluster = 'cluster_022_n_19'
working_dir = cluster+'_dir'
if not os.path.isdir(working_dir):
	os.mkdir(working_dir)

self_maps = 1

def map_and_merge(graph, fname1, fname2, out):
	os.system(f"minimap2 -x asm5 -D -c  {fname1} {fname2} > {out}.paf")
	paf = parse_paf(f'{out}.paf')
	paf.sort(key=lambda x:-x['aligned_bases'])
	merged_blocks = set()
	for hit in paf:
		if hit['query']['name'] in merged_blocks \
		   or hit['ref']['name'] in merged_blocks \
		   or hit['ref']['name']==hit['query']['name']:
			continue
		if set(graph.blocks[hit['query']['name']].sequences.keys()).intersection(graph.blocks[hit['ref']['name']].sequences.keys()):
			continue

		print('merging', fname1, fname2, hit)
		n.graph.merge_hit(hit)
		merged_blocks.add(hit['ref']['name'])
		merged_blocks.add(hit['query']['name'])


seqs = list(SeqIO.parse(f'inc_clusters/{cluster}.fasta','fasta'))
T = Phylo.read(f'inc_clusters/{cluster}.nwk','newick')

ni = 0
for seq in seqs:
	seq.id += f"_{ni:03d}"
	ni+=1

for (n,seq) in zip(T.get_terminals(), seqs):
	n.graph = Graph.from_sequence(seq.id, str(seq.seq).upper())
	n.name = seq.id
	n.fasta_fname = os.path.join(*[working_dir, n.name+'.fasta'])
	n.graph.to_fasta(n.fasta_fname)

node_count = 0
for n in T.get_nonterminals(order='postorder'):
	n.name = f'NODE_{node_count:07d}'
	node_count+=1
	n.graph = Graph.fuse(n.clades[0].graph, n.clades[1].graph)
	map_and_merge(n.graph, n.clades[0].fasta_fname, n.clades[1].fasta_fname, f"{working_dir}/{n.name}")
	n.fasta_fname = os.path.join(*[working_dir, n.name+'.fasta'])

	for i in range(self_maps):
		n.graph.to_fasta(n.fasta_fname+f'_iter_{i}')
		map_and_merge(n.graph, n.fasta_fname+f'_iter_{i}', n.fasta_fname+f'_iter_{i}', f"{working_dir}/{n.name}_iter_{i}")

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
				for j,o,r in zip(range(100), orig[i*100:(i+1)*100], rec[i*100:(i+1)*100]):
					if o!=r:
						print(i*100+j,o,r)
				node_path = T.get_path(n.name)
				for ni,pn in enumerate(node_path):
					print(ni, pn.name, orig==pn.graph.extract(n.name))
				import ipdb; ipdb.set_trace()
				print(i,orig[i*100:(i+1)*100])
				print(i,rec[i*100:(i+1)*100])
	print(len(T.root.graph.extract(n.name)))

T.root.graph.to_json(os.path.join(working_dir, 'graph.json'))
