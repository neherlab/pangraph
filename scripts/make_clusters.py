import numpy as np

from util import parse_kmer

# ------------------------------------------------------------------------
# Global constants/variables

graph_path = "/home/nolln/code/bio/pangraph/data/graphdist.nomap.npz"
kmers_path = "/home/nolln/code/bio/pangraph/data/kmerdist.txt"

# ------------------------------------------------------------------------
# Functions

def main():
    data = np.load(graph_path)
    graph_D, graph_names = data['arr_1'], data['arr_0']
    kmers_D, kmers_names = parse_kmer(kmers_path)

    matches = np.zeros(len(graph_names), dtype=np.int)
    for i, name in enumerate(graph_names):
        matches[i] = int(np.where(kmers_names == name)[0][0])

    tmp = np.zeros((len(graph_names), len(graph_names)))
    for i, mi in enumerate(matches):
        for j, mj in enumerate(matches):
            tmp[i, j] = kmers_D[mi, mj]
    kmers_D = tmp

    D = .5 * (graph_D + graph_D.T)
    D[kmers_D < .5] = np.inf

    return D

# ------------------------------------------------------------------------
# Main point of entry

if __name__ == "__main__":
    D = main()
