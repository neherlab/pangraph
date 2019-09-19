import argparse
import numpy as np

def to_newick(node, newick, parentdist, leaf_names):
    if node.is_leaf():
        return "%s:%.2f%s" % (leaf_names[node.id], parentdist - node.dist, newick)
    else:
        if len(newick) > 0:
            newick = "):%.2f%s" % (parentdist - node.dist, newick)
        else:
            newick = ");"
        newick = to_newick(node.get_left(), newick, node.dist, leaf_names)
        newick = to_newick(node.get_right(), ",%s" % (newick), node.dist, leaf_names)
        newick = "(%s" % (newick)

        return newick

def cluster_distance_matrix(fname):
    from scipy.spatial.distance import squareform
    from scipy.cluster.hierarchy import linkage, to_tree

    assert fname.endswith(".npz")

    data = np.load(fname)
    isos = data['arr_0']
    D    = data['arr_1']

    D[np.isinf(D)] = 200 #TODO: Fix this!
    D = .5 * (D + np.transpose(D))

    Dsqf = squareform(D)
    Z    = linkage(Dsqf, method="complete")
    T    = to_tree(Z)

    nwk = to_newick(T, "", T.dist, isos)
    with open(fname.replace(".npz", ".nwk"), "w") as fh:
        fh.write(nwk + "\n")

parser = argparse.ArgumentParser(description = "", usage = "map each isolate against entire database")
parser.add_argument("mtx",  type = str, help = "matrix file")

if __name__ == "__main__":
    args = parser.parse_args()
    cluster_distance_matrix(args.mtx)
