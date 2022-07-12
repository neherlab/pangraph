import numpy as np
import matplotlib.pyplot as plt

from scipy.optimize import curve_fit


def sharedness_block_mask(block_list, block_pa, present_in, absent_in):
    """
    given the genome to project onto, and the other strains of interest,
    returns a binary vector indicating whether block occurrences are shared or not by the strains of interest.
    """
    is_present = block_pa.loc[present_in].apply(np.all)
    is_absent = (~block_pa.loc[absent_in]).apply(np.all)
    is_accepted = is_present & is_absent
    block_mask = is_accepted[block_list]
    return block_mask


def block_length_vector(pan, strain):
    """
    Given the pangraph object and a particular strains, it returns a vector of
    block lengths, ordered as the block occurrences. Lengths keep into account
    blocks insertions and deletions
    """
    # extract path for specified strains
    path = pan.paths[strain]
    # list of block ids for the strain
    bl_ids = path.block_ids
    # initialize length container
    Ls = np.zeros(len(bl_ids), dtype=int)
    # build arrays for indexing block occurrences
    bl_strain, bl_num, bl_strand = (
        path.block_strains,
        path.block_nums,
        path.block_strands,
    )
    for n_idx, index in enumerate(zip(bl_strain, bl_num, bl_strand)):
        # capture block object
        bl_id = bl_ids[n_idx]
        bl = pan.blocks[bl_id]
        # capture block occurrence insertions and deletions
        ald = bl.align_diff[index]
        ins, delet = ald["insert"], ald["delete"]
        # count bp difference due to indels
        deltaL = 0
        for my_ins in ins:
            deltaL += len(my_ins[1])
        for my_del in delet:
            deltaL -= my_del[1]
        # insert corrected length
        Ls[n_idx] = len(bl) + deltaL

    return Ls


def produce_X(pan, main_strain, present_in, absent_in=[], block_pa_df=None):

    # if the block presence/absence dataframe is not provided then evaluate it
    if block_pa_df is None:
        block_pa_df = pan.to_blockcount_df() > 0

    block_ids = pan.paths[main_strain].block_ids
    block_mask = sharedness_block_mask(block_ids, block_pa_df, present_in, absent_in)
    Ls = block_length_vector(pan, main_strain)
    X = np.repeat(block_mask.values, Ls)

    return X


def exp_fit(x, p0, l):
    return p0 + (1.0 - p0) * np.exp(-x / l)


def check_autocorrelation(X, dx, xmax, visual_check=False):
    # fit with exponential + saturation

    L = X.size
    Ls = np.sum(X)
    # Lc = L - Ls

    V = np.arange(0, xmax, dx, dtype=int)
    Y = [np.sum(X * np.roll(X, v)) / Ls for v in V]

    pars, cov = curve_fit(exp_fit, V, Y, p0=(Y[-1], -2 * dx / (Y[2] - Y[0])))

    p0, l = pars[0], pars[1]

    if visual_check:
        fig, ax = plt.subplots(1, 1)
        axi = ax
        axi.plot(V, Y, "-")
        axi.plot(V, exp_fit(V, p0, l), "--")
        axi.set_xlabel("shift")
        axi.text(0.7, 0.9, f"p_0 = {p0:.3}", transform=axi.transAxes)
        axi.text(0.7, 0.85, f"l = {l:.3}", transform=axi.transAxes)

        # axi.text(0.7, 0.6, f"{str_i}", transform=axi.transAxes)
        # axi.text(0.7, 0.55, f"{str_j}", transform=axi.transAxes)

        axi.set_ylabel("Convolution")
        plt.tight_layout()
        plt.show()

    return V, Y, p0, l
