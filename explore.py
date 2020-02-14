import json
import pickle
import numpy as np

from glob import glob
from py.graph import Graph

cls   = "024"

betas = [1, 2, 5, 10, 25, 50, 75, 100, 250, 500, 750, 1000]
mus   = [1e0, 5e0, 1e1, 5e1, 1e2, 5e2, 1e3, 5e3, 1e4]

def getmu(f):
    nm = f.replace(f"data/graph/pkl/cluster_{cls}_m", "").split("_")[0]
    return int(nm), mus.index(int(nm))

def getbeta(f):
    nm = f.replace(f"data/graph/pkl/cluster_{cls}_m", "").split("_")[1]
    nm = nm.replace(".pkl", "")[1:]
    return int(nm), betas.index(int(nm))

if __name__ == "__main__":

    nblks = np.zeros((len(mus), len(betas)))
    xratio = np.zeros((len(mus), len(betas)))
    blklen = np.zeros((len(mus), len(betas)))

    Gs = []
    for f in glob(f"data/graph/pkl/cluster_{cls}*.pkl"):
        if f == f"data/graph/pkl/cluster_{cls}.pkl":
            continue

        mu,   imu = getmu(f)
        beta, ibt = getbeta(f)

        G = Graph.fromdict(pickle.load(open(f,'rb')))

        nblks[imu, ibt]  = len(G.blks)
        xratio[imu, ibt] = G.compressratio()
        blklen[imu, ibt] = sum(len(b.seq) for b in G.blks.values())/len(G.blks)

