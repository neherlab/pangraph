import numpy as np


class Branching:
    def __init__(self, pan, core_block_id):
        self.block_id = core_block_id
        self.blocks = []  # list of [bl, b, br] for every occurrence
        self.occs = []  # list of [ol, o, or]
        self.strains = []
        self.dir_r = []
        self.dir_l = []

        bl = pan.blocks[core_block_id]
        assert bl.frequency() == len(pan.paths), "must be a core block"
        assert bl.is_duplicated() == False, "must be a core block"

        for pt in pan.paths:
            blks, occs = find_next_blocks_from_path(pt, core_block_id)
            self.strains.append(pt.name)
            self.blocks.append(blks)
            self.occs.append(occs)
            same_r = occs[1][2] == occs[2][2]
            self.dir_r.append("_s" if same_r else "_d")
            same_l = occs[1][2] == occs[0][2]
            self.dir_l.append("_s" if same_l else "_d")

        self.blocks = np.array(self.blocks)
        self.dir_l = np.array(self.dir_l)
        self.dir_r = np.array(self.dir_r)

    def pick_side(self, side):
        """Return list of blocks, direction coincidence and occurrences
        for the corresponding side"""
        if side == "R":
            return (
                self.blocks[:, 2],
                self.dir_r,
                np.array([o[2] for o in self.occs]),
            )
        elif side == "L":
            return (
                self.blocks[:, 0],
                self.dir_l,
                np.array([o[0] for o in self.occs]),
            )
        else:
            raise Exception("side must be either `R` or `L`")

    def is_valid(self, side, pan, depth_lower_lim=2, len_lower_lim=300):
        """Function to determine whether the split is valid.
        The split is valid if two branches with the following characteristics are found:
        - outgoing depth >= `depth_lower_lim` (default 2)
        - length >= `len_lower_lim` (default 300 bp)"""

        bl_side, dir_side, occs = self.pick_side(side)
        bl_union = np.char.add(bl_side, dir_side)

        # count them, and arrange in decreasing count order
        res, ct = np.unique(bl_union, return_counts=True)
        srt = np.argsort(ct)[::-1]
        res, ct = res[srt], ct[srt]

        # select the two branches with highest count and satisfying the requirements
        selected = []
        for r, c in zip(res, ct):
            ct_condition = c >= depth_lower_lim
            L_condition = len(pan.blocks[r[:-2]]) >= len_lower_lim
            if ct_condition and L_condition:
                selected.append(r)
                if len(selected) > 1:
                    break

        if len(selected) < 2:
            return False, None

        results = []
        for bl_tag in selected:
            mask = bl_union == bl_tag
            results.append(
                {
                    "block": bl_tag[:-2],
                    "occs": occs[mask],
                    "direction": bl_tag[-2:],
                    "counts": mask.sum(),
                }
            )

        return True, results


def find_next_blocks_from_path(path, core_block_name):
    """given a path, finds the occurrence of the block and returns
    a list of three occurrences: the block ad its left and right neighbors"""
    L = len(path)
    w = np.argwhere(path.block_ids == core_block_name).flatten()
    assert len(w) == 1, "the block must be a core block"
    w = w[0]
    # previous and next block index
    b, n = (w - 1) % L, (w + 1) % L

    occs = []
    blocks = []
    for i in [b, w, n]:
        blocks.append(path.block_ids[i])
        occ = (path.name, path.block_nums[i], path.block_strands[i])
        occs.append(occ)

    # invert order if the middle block is on the other strand
    if not occs[1][2]:
        occs = occs[::-1]
        blocks = blocks[::-1]

    return blocks, occs
