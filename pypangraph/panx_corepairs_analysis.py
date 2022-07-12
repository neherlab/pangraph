# Functions that perform core-gene gaps analysis for panx data.
# The aim is to look at single-copy accessory genes difference
# in accessory genome between a specified pair of core genes that
# are always adjacent in all genomes.

import numpy as np
from collections import Counter


def coregaps_analysis(paths):
    """Perform coregaps analysis.
    Paths must be a dictionary {strain_id -> [list of gids]}.
    """

    # find core genes and duplicated genes
    core_genes = find_coregenes(paths)
    dupl_genes = find_duplicated_genes(paths)

    # raise error if too few core genes
    Ncg = len(core_genes)
    if Ncg <= 2:
        raise RuntimeError(f"The number of core-genes is too low: {Ncg}")

    # core-pairs gaps dictionary
    gap_dicts = {}
    for strain, p in paths.items():
        gap_dicts[strain] = build_corepair_dict(p, core_genes, dupl_genes)

    # take only common core-gene pairs
    core_pairs = [set(gap_d) for strain, gap_d in gap_dicts.items()]
    core_pairs = set.intersection(*core_pairs)

    # # raise error if too few pairs
    # Npairs = len(core_pairs)
    # if Npairs <= 2:
    #     raise RuntimeError(
    #         f'The number of core-gene pairs is too low: {Npairs}')
    #
    # # leave only omnipresent core-gene pairs
    # for str in gap_dicts:
    #     gap_dicts[str] = {cp: gap for cp, gap in gap_dicts[str].items()
    #                       if cp in core_pairs}

    return gap_dicts, core_pairs


def find_coregenes(paths):
    """Find single-copy core genes. Returns the results as a set."""
    sets = []
    for s, p in paths.items():
        # remove duplicated genes
        count = Counter(p)
        genes = {gid for gid, c in count.items() if c == 1}
        sets.append(genes)
    # core genes are the intersection of these sets.
    core_genes = set.intersection(*sets)
    return core_genes


def find_duplicated_genes(paths):
    """Find duplicated genes. Returns the results as a set."""
    dupl_genes = set()
    for s, p in paths.items():
        count = Counter(p)
        # take only genes that appear more than once
        d_genes = [gid for gid, c in count.items() if c > 1]
        dupl_genes.update(d_genes)
    return dupl_genes


def build_corepair_dict(path, core_genes, dupl_genes):
    """
    given a path and the set of coregenes, it builds a dictionary
    of core gene gaps and their accessory genes. Duplicated genes are not considered
    in the gap, and duplicated genes are removed.

    Args:
        - path (list) : ordered list of gene ids.
        - core_genes (set) : set of core genes.
        - dupl_genes (set) : set of duplicated genes.

    Returns:
        - gaps (dict) : dictionary core-pair -> accessory genes in the gap
            with form {(coregeneA, coregeneB) -> {accgeneA, accgeneB, ...}}

    Nb: the number of core-genes must be higher than 0 for the function to terminate,
    and higher than 2 to provide meaningful results
    """
    # core-pair gaps dictionary
    gaps = {}
    L, l = len(path), 0
    left_cg, right_cg = None, None
    first_cg = None  # first core gene found, used for exit condition
    gap_genes = set()  # accumulates non-core genes as we scroll the genome

    exit_condition = False
    while not exit_condition:
        gid = path[l]  # capture gene
        if gid in core_genes:  # if core gene
            if first_cg is None:  # if first, reset gap genes and set first_cg
                first_cg = gid
                gap_genes = set()
                left_cg = gid
            else:  # else add the gap and the core-gene pair
                right_cg = gid
                idx = sorted_tuple(left_cg, right_cg)
                gaps[idx] = gap_genes - dupl_genes  # remove duplicated genes
                gap_genes = set()  # reset
                left_cg = right_cg

                if right_cg == first_cg:  # if back to the first core-gene then exit
                    exit_condition = True
        else:  # if not core gene, extend the gap
            gap_genes.add(gid)
        l = (l + 1) % L
    return gaps


def sorted_tuple(A, B):
    if A <= B:
        return (A, B)
    else:
        return (B, A)


def pairwise_private_gapgenes(gaps_dict, core_pairs, strA, strB):
    """For a given pair of strains, it returns a list of differences
    in species-omnipresent core-gene gaps.

    Args:
        - gaps_dict (dict): dictionary {strain -> (core-gene pair) -> {accessory gap genes}}
            returned by the coregaps_analysis function
        - core_pairs (set): set of omni-species core-gene pairs, returned by
            the coregaps_analysis function
        - strA, strB (strings): specified strains

    Returns:
        - diff_gaps (dict): dictionary {(core-gene pair) -> {private genes}}
            for gaps with non-zero number of private genes.
    """

    diff_gaps = {}
    gapA, gapB = gaps_dict[strA], gaps_dict[strB]

    for cp in core_pairs:
        diff = gapA[cp] ^ gapB[cp]  # symmetric difference
        if len(diff) > 0:
            diff_gaps[cp] = diff

    return diff_gaps
