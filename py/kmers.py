import os, sys
import numpy as np
import json

from Bio import Phylo, Seq, SeqIO, SeqRecord
from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor

from utils import Strand, asstring, parsepaf, panic, tryprint, asrecord, newstrand
from graph import Graph

# ------------------------------------------------------------------------
# Global variables

maxselfmaps = 25

# ------------------------------------------------------------------------
# Helper functions

def flatten(x):
    return np.ndarray.flatten(x)

def nodiag(mtx):
    return mtx-np.diag(np.diag(mtx))

def parse(mtx):
    with open(mtx) as fh:
        nrows = int(fh.readline().strip())
        M, r  = np.zeros((nrows, nrows), dtype=float), 0

        del_idxs  = []
        seq_names = []
        for li, line in enumerate(fh):
            e = line.strip().split()
            n = e[0].split('/')[-1][:-3]
            if n not in seq_names:
                seq_names.append(n)
                M[li,:(li+1)] = [float(x) for x in e[1:]]
            else:
                del_idxs.append(li)

    M = np.delete(M, del_idxs, axis=0)
    M = np.delete(M, del_idxs, axis=1)

    # Symmetrize
    M = 1 - M
    M = nodiag(M + M.T)/2

    return M, seq_names

def tolist(dmtx):
    assert len(dmtx.shape) == 2 and dmtx.shape[0] == dmtx.shape[1], "expected a square matrix"

    dlst = []
    for n in range(dmtx.shape[0]):
        dlst.append(list(dmtx[n,:(n+1)]))

    return dlst

# ------------------------------------------------------------------------
# Node and Tree classes

class Node(object):

    # ---------------------------------
    # Internal functions

    def __init__(self, name, parent, dist, children=[]):
        self.name     = name
        self.dist     = dist
        self.parent   = parent
        self.children = children
        self.fapath   = ""
        self.graph    = None

    def __str__(self):
        if self.dist is None:
            return f"{self.name} :: Unknown"
        else:
            return f"{self.name} :: {self.dist:.4f}"

    def __repr__(self):
        return self.__str__()

    # ---------------------------------
    # Static functions

    @classmethod
    def fromdict(cls, d, parent):
        N = Node(d['name'], parent, d['dist'])
        N.children = [Node.fromdict(child, N) for child in d['children']]
        N.fapath = d['fapath']
        N.graph  = Graph.fromdict(d['graph']) if d['graph'] is not None else None

        return N

    # ---------------------------------
    # Class methods

    def isleaf(self):
        return len(self.children) == 0

    def postorder(self):
        for child in self.children:
            for it in child.postorder():
                yield it
        yield self

    def newparent(self, parent, dist):
        self.parent = parent
        self.dist   = dist if dist > 0 else 0

    def tonwk(self, wtr):
        if not self.isleaf():
            wtr.write("(")
            for i, child in enumerate(self.children):
                if i > 0:
                    wtr.write(",")
                child.tonwk(wtr)
            wtr.write(")")

        wtr.write(self.name)
        wtr.write(":")
        wtr.write(f"{self.dist:.6f}")

    def tojson(self):
        return {'name'     : self.name,
                'dist'     : self.dist,
                'children' : [ child.tojson() for child in self.children ],
                'fapath'   : self.fapath,
                'graph'    : self.graph.todict() if self.graph is not None else None }

class Tree(object):

    # ------------------- 
    # Class constructor

    def __init__(self, bare=False):
        self.root   = Node("ROOT", None, 0) if not bare else None
        self.seqs   = None
        self.leaves = None

    # ------------------- 
    # Static methods

    # Loading from json
    @classmethod
    def fromjson(cls, rdr):
        data   = json.load(rdr)
        T      = Tree(bare=True)
        T.root = Node.fromdict(data['tree'], None)

        return T

    # Neighbor joining
    # BioPython implementation is too slow.
    @classmethod
    def nj(cls, mtx, names, verbose=False):
        assert len(names) == len(set(names)), "non-unique names found"

        T = Tree()
        for name in names:
            T.root.children.append(Node(name, T.root, None, children=[]))
        idx = 0

        # -- Internal functions --
        def calcq(D):
            n = D.shape[0]
            Q = (n-2)*D - np.sum(D,axis=1) - np.sum(D,axis=0)
            np.fill_diagonal(Q, np.inf)

            return Q

        def minpair(q):
            i, j = np.unravel_index(np.argmin(q), q.shape)
            qmin = q[i,j]
            if i > j:
                i, j = j, i
            return (i,j), qmin

        def pairdists(D, i, j):
            n  = D.shape[0]
            d1 = .5*D[i,j] + 1/(2*(n-2)) * (np.sum(D[i,:], axis=0) - np.sum(D[j,:], axis=0))
            d2 = D[i,j] - d1
            return d1, d2

        def newdists(D, i, j):
            return .5*(D[i, :] + D[j,:] - D[i,j])

        def join(D, debug=False):
            nonlocal idx
            Q = calcq(D)
            (i, j), qmin = minpair(Q)
            assert i < j
            if debug:
                q0min = min(flatten(Q[:]))
                assert abs(qmin-q0min) < 1e-2, f"minimum not found correctly. returned {qmin}, expected {q0min}"
                print(f"{D}\n--> Joining {i} and {j}. d={D[i,j]}")

            node   = Node(f"NODE_{idx:05d}", T.root, None, [T.root.children[i], T.root.children[j]])

            d1, d2 = pairdists(D, i, j)
            dnew   = newdists(D, i, j)
            node.children[0].newparent(node, d1)
            node.children[1].newparent(node, d2)

            D[i, :] = dnew
            D[:, i] = dnew
            D[i, i] = 0
            D = np.delete(D, j, axis=0)
            D = np.delete(D, j, axis=1)
            T.root.children[i] = node
            T.root.children.pop(j)

            idx = idx + 1

            return D

        while mtx.shape[0] > 2:
            if verbose:
                print(f"--> Matrix size={mtx.shape[0]}. Number of root children={len(T.root.children)}")
            mtx = join(mtx)

        assert mtx.shape[0] == 2
        d = mtx[0, 1]
        T.root.children[0].dist = d/2
        T.root.children[1].dist = d/2

        return T

    # ------------------- 
    # Class methods 

    def postorder(self):
        return self.root.postorder()

    def getleafs(self):
        if self.leaves is None:
            self.leaves = [node for node in self.postorder() if node.isleaf()]

        return self.leaves

    def numleafs(self):
        if self.leaves is None:
            self.leaves = [node for node in self.postorder() if node.isleaf()]

        return len(self.leaves)

    def align(self, seqs, save=True, verbose=False):
        # Debugging function that will check reconstructed sequence against known real one.
        def check(seqs, T, G, verbose=False):
            nerror = 0
            uncompressed_length = 0
            for n in T.getleafs():
                if n.name not in G.seqs:
                    continue

                seq  = seqs[n.name]
                orig = str(seq[:].seq).upper()
                tryprint(f"--> Checking {n.name}", verbose=verbose)
                rec  = G.extract(n.name)
                uncompressed_length += len(orig)
                if orig != rec:
                    nerror += 1

                    with open("test.fa", "w+") as out:
                        out.write(f">original\n{orig}\n")
                        out.write(f">reconstructed\n{rec}")

                    for i in range(len(orig)//100):
                        if (orig[i*100:(i+1)*100] != rec[i*100:(i+1)*100]):
                            print("-----------------")
                            print("O:", i, orig[i*100:(i+1)*100])
                            print("G:", i, rec[i*100:(i+1)*100])

                            diffs = [i for i in range(len(rec)) if rec[i] != orig[i]]
                            pos   = [0]
                            blks  = G.seqs[n.name]
                            for b, strand, num in blks:
                                pos.append(pos[-1] + len(G.blks[b].extract(n.name, num)))
                            pos = pos[1:]

                            testseqs = []
                            for b in G.seqs[n.name]:
                                if b[1] == Strand.Plus:
                                    testseqs.append("".join(G.blks[b[0]].extract(n.name, b[2])))
                                else:
                                    testseqs.append("".join(Seq.reverse_complement(G.blks[b[0]].extract(n.name, b[2]))))

                else:
                    tryprint(f"+++ Verified {n.name}", verbose=verbose)

            if nerror == 0:
                tryprint("all sequences correctly reconstructed", verbose=verbose)
                tlength = np.sum([len(x) for x in G.blks.values()])
                tryprint(f"--- total graph length: {tlength}", verbose=verbose)
                tryprint(f"--- total input sequence: {uncompressed_length}", verbose=verbose)
                tryprint(f"--- compression: {uncompressed_length/tlength:1.2f}", verbose=verbose)
            else:
                raise ValueError("bad sequence reconstruction")

        # T = Phylo.read(path, "newick")
        if self.numleafs() == 1:
            return Graph()

        for i, n in enumerate(self.getleafs()):
            seq      = asrecord(seqs[n.name])
            n.graph  = Graph.fromseq(seq.id, str(seq.seq).upper())
            n.name   = seq.id
            n.fapath = f"{Graph.blddir}/{n.name}"
            tryprint(f"------> Outputting {n.fapath}", verbose=verbose)
            n.graph.tofasta(n.fapath)

        nnodes = 0
        for n in self.postorder():
            try:
                if n.isleaf():
                    continue
                nnodes += 1

                # Simple graph "fuse". Straight concatenation
                print(f"Fusing {n.children[0].name} with {n.children[1].name}", file=sys.stderr)
                n.graph  = Graph.fuse(n.children[0].graph, n.children[1].graph)
                # check(seqs, T, n.graph)
                n.fapath = os.path.join(*[Graph.blddir, f"{n.name}.fasta"])

                tryprint(f"-- node {n.name} with {len(n.children)} children", verbose)

                # Initial map from graph to itself
                n.graph, _ = n.graph._mapandmerge(n.children[0].fapath, n.children[1].fapath,
                                        f"{Graph.blddir}/{n.name}")

                i, contin = 0, True
                while contin:
                    tryprint(f"----> merge round {i}", verbose)
                    # check(seqs, T, n.graph)
                    itrseq = f"{Graph.blddir}/{n.name}_iter_{i}"
                    n.graph.tofasta(itrseq)
                    n.graph, contin = n.graph._mapandmerge(itrseq, itrseq, f"{Graph.blddir}/{n.name}_iter_{i}")
                    i += 1

                    contin &= i < maxselfmaps

                tryprint(f"-- Blocks: {len(n.graph.blks)}, length: {np.sum([len(b) for b in n.graph.blks.values()])}\n", verbose)
                n.graph.tofasta(f"{Graph.blddir}/{n.name}")
                # Continuous error logging
                print(f"Node {n.name}", file=sys.stderr)
                print(f"--> Compression ratio parent: {n.graph.compressratio()}", file=sys.stderr)
                print(f"--> Compression ratio child1: {n.children[0].graph.compressratio()}", file=sys.stderr)
                print(f"--> Compression ratio child2: {n.children[1].graph.compressratio()}", file=sys.stderr)
                # Output with content
                print(f"{n.graph.compressratio()}\t{n.children[0].graph.compressratio()}\t{n.children[1].graph.compressratio()}\t{n.children[0].dist+ n.children[1].dist}", file=sys.stdout, flush=True)
            except:
                print(f"ERROR AT NODE {n.name}")
                with open("error.json", "w+") as out:
                    self.tojson(out)

    def tonwk(self, wtr):
        self.root.tonwk(wtr)
        wtr.write(";")

    def tojson(self, wtr):
        data = {'tree' : self.root.tojson()}
        wtr.write(json.dumps(data))

# ------------------------------------------------------------------------
# Main point of entry for testing

if __name__ == "__main__":
    # M   = np.array([[0, 2, 4, 4], [2, 0, 4, 4], [4, 4, 0, 2], [4, 4, 2, 0.0]])
    # nms = ["A", "B", "C", "D"]
    M, nms = parse("data/kmerdist.txt")
    T = Tree.nj(M, nms)

    import pyfaidx as fai
    seqs = fai.Fasta("data/all_plasmids_filtered.fa")
    T.align(seqs)

    print("DUMPING")
    with open("data/tree.json", "w+") as fd:
        T.tojson(fd)

    S1 = json.dumps(T.root.tojson())
    print("LOADING")
    with open("data/tree.json", "r") as fd:
        T = Tree.fromjson(fd)

    S2 = json.dumps(T.root.tojson())

    print(S1==S2)
