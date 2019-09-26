import numpy as np
import glob
from Bio import SeqIO
from scipy.cluster.hierarchy import dendrogram, linkage, to_tree
from scipy.spatial.distance import squareform

def getNewick(node, newick, parentdist, leaf_names):
    if node.is_leaf():
        return "%s:%.2f%s" % (leaf_names[node.id], parentdist - node.dist, newick)
    else:
        if len(newick) > 0:
            newick = "):%.2f%s" % (parentdist - node.dist, newick)
        else:
            newick = ");"
        newick = getNewick(node.get_left(), newick, node.dist, leaf_names)
        newick = getNewick(node.get_right(), ",%s" % (newick), node.dist, leaf_names)
        newick = "(%s" % (newick)
        return newick

with open('kmer_dist_matrix.dist') as fh:
	nrows = int(fh.readline().strip())
	M = np.zeros((nrows, nrows), dtype=float)
	seq_names = []
	for li,line in enumerate(fh):
		e = line.strip().split()
		seq_names.append(e[0].split('/')[-1][:-3])
		M[li,:(li+1)] = [float(x) for x in e[1:]]

	Msym = M + M.T
	Msym -= np.eye(nrows)

distance = 1-Msym

cluster_list = glob.glob('inc_clusters/*fasta')

for c in cluster_list:
	cluster_seq_names = [s.id for s in SeqIO.parse(c, 'fasta')]
	unique_seqs = {s for s in cluster_seq_names}
	if len(unique_seqs)>1:
		indices = np.array(sorted([seq_names.index(x) for x in unique_seqs]))
		sub_matrix = squareform(distance[indices][:,indices])
		Z = linkage(sub_matrix)
		T = to_tree(Z)
		nwk = getNewick(T, "", T.dist, [seq_names[i] for i in indices])
		with open(c[:-6]+'.nwk', 'w') as fh:
			fh.write(nwk+'\n')
