# Functions to locate positions or intervals in the genome
# on blocks in pangraph

import numpy as np


class Locator:
    """Given a pangraph, builds a map that can be used to quickly
    locate a genomic position on the pangraph.
    """

    def __init__(self, pan):
        # build a map
        self.map = build_map(pan.paths, pan.blocks)

    def find_position(self, strain, pos):
        """Returns the block-id associated to a particular position
        on the genome. Position must be provided with 1-based indexing.
        It also returns the position of the nucleotide in the block sequence
        (1-based indexing), along with the block occurrence tag
        (strain, block n.,strand).
        Nb: the position is relative to the block sequence in the strain
        considered, and not to the block consensus.
        """
        pmap = self.map[strain]
        return pmap.position_to_block(pos)

    def find_interval(self, strain, pos_b, pos_e):
        """Returns the list of blocks corresponding to an interval
        of positions on the genomes.
        Positions should be provided with 1-based indexing.
        """
        pmap = self.map[strain]
        return pmap.interval_to_blocks(pos_b, pos_e)

    def __getitem__(self, idx):
        return self.map[idx]

    def __iter__(self):
        return iter(self.map)

    def items(self):
        return self.map.items()


class PathMap:
    def __init__(self, bl_ids, bl_begs, bl_ends, bl_lengths, bl_occs):
        """Takes care of reordering the lists of blocks based on beginning
        positions on the genomes.
        """
        self.ids = np.array(bl_ids)
        self.b = np.array(bl_begs)
        self.e = np.array(bl_ends)
        self.N = len(self.ids)
        self.Ls = np.array(bl_lengths)
        self.occs = np.array(bl_occs, dtype=object)
        self.path_L = np.sum(bl_lengths)
        assert np.all((self.e + 1 - self.b) % self.path_L == self.Ls)  # TODO: remove

        # reorder so that list of beginnings of blocks is ordered
        # computationally intensive
        order = np.argsort(self.b)
        self.b = self.b[order]
        self.e = self.e[order]
        self.ids = self.ids[order]
        self.Ls = self.Ls[order]
        self.occs = self.occs[order]

    def position_to_block_idx(self, pos):
        """Given a position on the genome, returns the index of the block
        in the PathMap lists.
        """
        idx = np.searchsorted(self.b, pos, side="right")
        idx = (idx - 1) % self.N  # correct for periodic boundary conditions
        return idx

    def position_to_block(self, pos):
        """Relates a position on the genome (1-based indexing!) to a position in
        a block. It returns the block id, the position of the nucleotide in the block
        (1-based indexing!) and the block occurrence tuple: (strain, block n., strand).
        """
        idx = self.position_to_block_idx(pos)
        bl_id = self.ids[idx]
        b, e, s, pthL = self.b[idx], self.e[idx], self.occs[idx][2], self.path_L
        bl_pos = position_in_block_coordinates(pos, b, e, s, pthL)
        occ = tuple(self.occs[idx])
        return bl_id, bl_pos, occ

    def interval_to_blocks(self, pos_b, pos_e):
        """Given a beginning and end positions (1-based indexing) on the genome,
        this function returns the pangraph blocks that contain this interval.
        For each block, it also returns the position of the gene relative
        to the block (1-based indexing). This position is relative to the
        particular block sequence, not the consensus. It also returns the corresponding
        list of block occurrences.
        """
        # find indices of start and end blocks
        idx_b = self.position_to_block_idx(pos_b)
        idx_e = self.position_to_block_idx(pos_e)

        # find interval start and end position in block frame of reference
        _, pb, _ = self.position_to_block(pos_b)
        _, pe, _ = self.position_to_block(pos_e)

        # create the list of indices of blocks that contain the interval
        wrap_1 = idx_e < idx_b
        wrap_2 = idx_e == idx_b
        strand = self.occs[idx_b][2]
        wrap_2 &= (strand & (pe < pb)) | ((not strand) & (pe > pb))
        if wrap_2:
            message = "warning: interval starts and ends in the same block"
            message += " but it wraps around the genome.\n"
            message += f"strand = {strand}, beg = {pos_b}, end = {pos_e},\n"
            message += f"block beg = {self.ids[idx_b]}, block end = {self.ids[idx_e]}"
            message += f"beg pos in block = {pb}, end pos in block = {pe}."
            print(message)
        if wrap_1 | wrap_2:
            idxs = np.arange(idx_b, idx_e + self.N + 1) % self.N
        else:
            idxs = np.arange(idx_b, idx_e + 1)

        I = [(1, self.Ls[idx]) for idx in idxs]
        bl_ids = self.ids[idxs]
        occs = self.occs[idxs]
        occs = [tuple(occ) for occ in occs]

        # set beginning and end
        Ib, occb = I[0], occs[0]
        I[0] = (pb, Ib[1]) if occb[2] else (Ib[0], pb)
        Ie, occe = I[-1], occs[-1]
        I[-1] = (Ie[0], pe) if occe[2] else (pe, Ie[1])
        return bl_ids, I, occs


def build_map(paths, blocks):
    """
    Given a set of paths and blocks, builds a dictionary {strain : PathMap} for every strain.
    """
    pan_map = {}
    for path in paths:
        strain = path.name
        bl_ids, bl_starts, bl_ends, bl_Ls, bl_occs = [[] for _ in range(5)]
        for bl, n, st in zip(path.block_ids, path.block_nums, path.block_strands):
            # index for the alignment dictionary in block object
            occ = (strain, n, st)

            # nb: julia indexing
            aln = blocks[bl].alignment
            pos_b, pos_e = aln.pos[occ]
            bl_L = aln.block_occurrence_length(occ)

            # build lists of block indices and positions
            bl_ids.append(bl)
            bl_starts.append(pos_b)
            bl_ends.append(pos_e)
            bl_Ls.append(bl_L)
            bl_occs.append(occ)

        # use the lists to build a map object for every strain
        pan_map[strain] = PathMap(bl_ids, bl_starts, bl_ends, bl_Ls, bl_occs)

    return pan_map


def position_in_block_coordinates(pos, bl_b, bl_e, bl_s, pth_L):
    """Given a position in the genome, information on the block hosting
    this position, returns the (1-based) index of the position in the block.
    The other arguments are:
    - block beginning/end position in the genome.
    - strandedness of block occurrence (forward/reverse strand)
    - total genome length
    All positions are in 1-based indexing.
    """
    assert (pos - bl_b) % pth_L <= (bl_e - bl_b) % pth_L
    assert (bl_e - pos) % pth_L <= (bl_e - bl_b) % pth_L
    if bl_s:
        return (pos - bl_b + 1) % pth_L
    else:
        return (bl_e - pos + 1) % pth_L
