import matplotlib.pyplot as plt
import numpy as np
from collections import Counter, defaultdict

from . import pangraph_interface as pgi

# TODO: implement distance measure in pairwise projection


class PanProjector:
    """Class to create and store pairwise projections of the
    pangraph on two particular strains.
    """

    def __init__(self, pan: pgi.Pangraph):

        # store a reference to the pangraph
        self.pan = pan

        # path dictionary
        self.paths = pan.to_paths_dict()

        # block statistics dataframe
        self.bls = pan.to_blockstats_df()

        # block lengths dictionary
        self.lengths = self.bls["len"].to_dict()

        # dictionary of duplicated genes
        N_dupl = self.bls["duplicated"].sum()
        self.is_dupl = self.bls["duplicated"].to_dict()

        # assign random color to duplicated genes
        dupl_mask = self.bls["duplicated"]
        cmap = plt.get_cmap("jet")
        dupl_blocks = self.bls.index[dupl_mask]
        self.dupl_color = {
            bl: cmap(np.random.rand()) for n, bl in enumerate(dupl_blocks)
        }

        # dictionary to store pairwise projections
        self.projections = {}

    def project(self, strA, strB, store=False, exclude_dupl=False):
        """Given a pair of strains, it builds the pairwise projection

        Args:
            - strA, strB (strings): name of the strands
        """

        pthA, pthB = self.pan.paths[strA], self.pan.paths[strB]

        # ordered list of blocks
        blA, blB = pthA.block_ids, pthB.block_ids

        # list of strandedness
        sA, sB = pthA.block_strands, pthB.block_strands

        # list of block lengths in nucleotides
        LsA = [self.lengths[b] for b in blA]
        LsB = [self.lengths[b] for b in blB]

        # exclude blocks that are duplicated in any genome
        if exclude_dupl:
            dupl = set(self.bls.index[self.bls["duplicated"]])
        else:
            dupl = None

        # evaluate projected graph
        proj_gr = ProjectedGraph(
            strA, strB, blA, blB, sA, sB, LsA, LsB, exclude_dupl=dupl
        )

        if store:
            # store projection into object's dictionary
            label = (strA, strB)
            self.projections[label] = proj_gr

        return proj_gr


def stitch_paths(blA, blB, strandA, strandB, exclude_dupl=None):
    """Function that creates a list of maximally extended
    agreeing parts of the two paths.
    """

    # capture lists of blocks and directions
    bA, bB = list(blA), list(blB)
    sA, sB = list(strandA), list(strandB)

    # list of blocks that are already taken
    is_taken_A = np.zeros(len(blA), dtype=bool)
    is_taken_B = np.zeros(len(blB), dtype=bool)

    # list of which duplicate is the same
    dupl_id_A = np.zeros(len(blA), dtype=int)
    dupl_id_B = np.zeros(len(blB), dtype=int)
    dupl_count = defaultdict(lambda: 1)

    # set of blocks that are duplicated either in A or B
    CA, CB = Counter(blA), Counter(blB)
    duplA = {x for x in CA if CA[x] > 1}
    duplB = {x for x in CB if CB[x] > 1}
    dupl = duplA | duplB
    if exclude_dupl is not None:
        dupl = dupl | exclude_dupl

    # set of seed blocks: common to both and not duplicated in either
    seeds = (set(blA) & set(blB)) - dupl

    # dictionary of results:
    res = defaultdict(list)

    while len(seeds) > 0:
        # wile seed blocks are still present, grow an extendable block from them.
        seed = seeds.pop()
        # retrieve seed index and strandedness
        i0A, i0B = bA.index(seed), bB.index(seed)
        sA0, sB0 = sA[i0A], sB[i0B]

        # create metablock from seed
        EB = ExtendableBlock(seed, sA0, sB0, i0A, i0B, len(bA), len(bB))

        # set starting block as taken
        is_taken_A[i0A], is_taken_B[i0B] = True, True

        while not EB.is_over:
            idxA, idxB = EB.next_candidate_idxs()
            cAi, cBi = bA[idxA], bB[idxB]
            sAi, sBi = sA[idxA], sB[idxB]
            takA, takB = is_taken_A[idxA], is_taken_B[idxB]
            success = EB.try_candidate(cAi, cBi, sAi, sBi, takA, takB)
            if success:
                is_taken_A[idxA], is_taken_B[idxB] = True, True
                if cAi in dupl:
                    dc = dupl_count[cAi]
                    dupl_count[cAi] += 1
                    dupl_id_A[idxA] = dc
                    dupl_id_B[idxB] = dc

        # append once growth is finished
        EB.append_results(res)
        # remove common blocks from seed set
        seeds = seeds - set(EB.seq)

    return res, is_taken_A, is_taken_B, dupl_id_A, dupl_id_B


class ExtendableBlock:
    """Extendable block object: starting from a seed
    point, it can be extended in both directions as long as the
    two paths agree.
    """

    def __init__(self, seed_bl, strandA, strandB, i0A, i0B, LtotA, LtotB):
        self.seq = [seed_bl]
        self.same = strandA == strandB
        self.stA, self.stB = [strandA], [strandB]
        self.i0A, self.i0B = i0A, i0B
        self.L = 1
        self.LtotA, self.LtotB = LtotA, LtotB
        self.is_over = False
        self.right = True

    def next_candidate_idxs(self):
        iA, iB = None, None
        if self.right:
            iA = self.i0A + self.L
        else:
            iA = self.i0A - 1

        if self.right ^ self.same:
            iB = self.i0B - 1
        else:
            iB = self.i0B + self.L

        return iA % self.LtotA, iB % self.LtotB

    def extend(self, block, strandA, strandB):
        self.L += 1
        self.stA.append(strandA)
        self.stB.append(strandB)

        if self.right:
            self.seq.append(block)
        else:
            self.seq.insert(0, block)
            self.i0A = (self.i0A - 1) % self.LtotA

        if self.right ^ self.same:
            self.i0B = (self.i0B - 1) % self.LtotB

    def try_candidate(self, cA, cB, strA, strB, takA, takB):

        # if one of the two is already taken, then failure
        if takA or takB:
            self.failed_direction()
            return False

        success = False
        if cA == cB:
            if self.same:
                if strA == strB:
                    success = True
            else:
                if strA != strB:
                    success = True

        if not success:
            self.failed_direction()
        else:
            self.extend(cA, strA, strB)

        return success

    def failed_direction(self):
        if self.right:
            self.right = False
        else:
            self.is_over = True

    def append_results(self, res):
        res["L"].append(self.L)
        res["i0A"].append(self.i0A)
        res["i0B"].append(self.i0B)
        res["seq"].append(self.seq)
        res["strA"].append(self.stA)
        res["strB"].append(self.stB)
        res["same_strand"].append(self.same)


class ProjectedGraph:
    """Class to create a graph projection"""

    def __init__(self, strA, strB, pthA, pthB, sA, sB, LsA, LsB, exclude_dupl):

        self.strA = strA
        self.strB = strB

        # evaluate common parts of the path
        st, commA, commB, did_A, did_B = stitch_paths(pthA, pthB, sA, sB, exclude_dupl)

        # define the two meta-paths
        self.MPA = MetaPath(strA, pthA, sA, LsA, did_A)
        self.MPB = MetaPath(strB, pthB, sB, LsB, did_B)

        # add common parts
        self.MPA.add_stitches(st, commA, "A")
        self.MPB.add_stitches(st, commB, "B")

        # finalize path creation
        self.MPA.metapath_creation()
        self.MPB.metapath_creation()


class MetaPath:
    """Object containing path and projected path representations"""

    def __init__(self, strain, pth, strand, Ls, did):
        self.pth = list(pth)
        self.L = len(pth)
        self.s = strand
        self.comm = np.zeros(self.L, dtype=bool)
        self.chunk_id = np.zeros(self.L, dtype=int)
        self.bl_Ls = np.array(Ls)
        self.dupl_id = did

    def add_stitches(self, st, comm, strain):

        idx0 = st[f"i0{strain}"]
        Ls = st["L"]
        seqs = st["seq"]

        K = len(Ls)

        for k in range(K):
            beg, end = idx0[k], idx0[k] + Ls[k]
            beg, end = beg % self.L, end % self.L

            if beg < end:
                self.comm[beg:end] = True
                self.chunk_id[beg:end] = k + 1
            else:
                self.comm[beg:], self.comm[:end] = True, True
                self.chunk_id[beg:], self.chunk_id[:end] = k + 1, k + 1

        assert np.all(self.comm == np.array(comm))

    def metapath_creation(self):
        # finalizing steps. Add negative chunk number to non-common chunks
        if np.sum(~self.comm) == 0:
            return

        current_id = -1
        current_chunk = []
        terminate = False
        l = np.argmax(self.comm)
        first_idx = None
        while not terminate:
            # print(l, current_id, first_idx, current_chunk)
            if self.comm[l] == 0:
                current_chunk.append(l)
                if first_idx is None:
                    first_idx = l
                elif first_idx == l:
                    terminate = True
            elif len(current_chunk) > 0:
                self.chunk_id[current_chunk] = current_id
                current_id -= 1
                current_chunk = []

            l = (l + 1) % self.L
