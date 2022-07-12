import numpy as np
import pandas as pd


class blockpairs_analysis:
    """Class that hanalizes the pangraph in terms of
    single-copy-core-block pairs, and the accessory blocks
    in between.
    """

    def __init__(self, pan):
        """
        Attributes of the class are:
        - pan (Pangraph object)
        - df (dataframe) : keys are block ids and columns are block
            properties
        - cgp (dict) : nested dicionary containing information on the core-block
            gaps for each strain. It has the form:
            {strain -> (bl_A, bl_B) -> set([interblocks])}
            Where core-block pairs are odered alphabetically
        - core_pairs (list) : list of core-block pairs that occurr in every strain.
            They can be used as index for the second layer of cgp.
        """

        # capture pangraph object
        self.pan = pan

        # creates dataframe of block properties
        # this df can be used to retrieve block lengths,
        # whether they are core, and whether they are duplicated
        self.df = pan.to_blockstats_df()

        # add frequency property
        S = len(pan.paths)  # number of strains
        self.df["freq"] = self.df["n. strains"] / S  # strain frequency

        # # create path dictionary and list of strains
        # self.paths = pan.to_paths_dict()
        # self.strains = list(self.paths.keys())

        # generates self.cgp and self.core_pairs
        self.build_coregenepair_dict()

    def build_coregenepair_dict(self):
        """Function that builds a nested dictionary `self.cgp`
        {strain -> (coreblock_A, coreblock_B) -> set([bl1, bl2 ...])}
        where the list contains blocks that occurr between the
        core-block pairs.
        It also creates a list `self.core_pairs` containing the
        core-block pairs that are common to all individuals
        """

        # initialize the dictionary
        self.cgp = {}

        # create dictionary that links each block to whether it is a core block
        is_coreblock = self.df["core"].to_dict()

        # raise error if too few core blocks for this analysis
        tot_n_coreblocks = self.df["core"].sum()
        if tot_n_coreblocks <= 2:
            raise ValueError(
                f"Error: only {tot_n_coreblocks} core blocks in the graph.\
                                not enough for a pairwise analysis."
            )

        # build core block pairs dictionary
        paths = self.pan.to_paths_dict()
        for strain, path in paths.items():
            self.cgp[strain] = parse_coreblock_pairs(path, is_coreblock)

        # generate a list of common pairs
        self.core_pairs = sieve_commonkeys(list(self.cgp.values()))

    def private_block_list(self, strain_A, strain_B):
        """Given a pair of strains, it returns a dictionary
        {pair -> set([interblocks])} listing private blocks for every
        core-block pair. Private blocks are defined as blocks that are
        present in one strain but not the other.
        """

        gaps = {}  # list of interblocks sets. One set per coregene-pair gap
        for pair in self.core_pairs:
            A = self.cgp[strain_A][pair]
            B = self.cgp[strain_B][pair]
            gaps[pair] = A ^ B  # block that are only in either A or B

        return gaps


def parse_coreblock_pairs(path, is_core):
    """Given a path and a dictionary {bl : is_core}, this function returns
    a dictionary of the form {(corebl_A, corebl_B) : set([inter-blocks])}
    NB: there must be more than 2 core genes in the path, otherwise
    results are meaningless.
    """
    Lbl, Rbl = None, None  # containers for left and right coreblocks
    firstblock = None  # save first coreblock, to check when we are
    P = len(path)  # total number of blocks in the path
    error_flag = (
        False  # flag to raise an error in case the path is parsed more than twice
    )

    corepair_dict = {}  # initialize the dictionary
    interblocks = []  # container for inter-blocks
    i = 0  # path index
    while (Rbl is None) or (
        Rbl != firstblock
    ):  # stop when it loops back to the first core block
        bl = path[i]
        if is_core[bl]:
            if Lbl is None:  # first core block found
                interblocks = []
                Lbl = bl
                firstblock = bl
            else:  # non-first core block found
                Rbl = bl
                idx = sorted_pair(Lbl, Rbl)
                # save as set (removes duplicated blocks!)
                corepair_dict[idx] = set(interblocks)
                interblocks = []
                Lbl = Rbl
        else:
            interblocks.append(bl)

        # check how many times periodic boundary conditions are passed
        if i + 1 == P:
            if error_flag:
                raise ValueError(
                    "The path was already parsed twice. Probably too few\
                                 core blocks?"
                )
            else:
                error_flag = True

        # periodic boundary conditions
        i = (i + 1) % P

    return corepair_dict


def sorted_pair(A, B):
    if A <= B:
        return (A, B)
    else:
        return (B, A)


def sieve_commonkeys(dict_list):
    """given a list of dictionaries, this function returns
    the list of all keys that are common to all dictionaries.
    """
    assert isinstance(dict_list, list), "dictionaries must be passed as a list"

    common_keys = set(dict_list[0].keys())

    for dct in dict_list[1:]:
        common_keys &= set(dct.keys())

    return list(common_keys)
