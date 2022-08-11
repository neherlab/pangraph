from collections import defaultdict
import json
import glob
import re
import numpy as np
import itertools as itt

# list of species
species = [
    "klebsiella_pneumoniae",
    "helicobacter_pylori",
    "prochlorococcus_marinus",
    "mycobacterium_tuberculosis",
    "escherichia_coli",
]

# list of accession numbers to exclude for each species
exclude = {
    "klebsiella_pneumoniae": ["NZ_CP012300", "NC_011283", "NZ_AP014950", "NZ_CP016811"],
    "helicobacter_pylori": [
        "NC_022130",
        "NC_017357",
        "NC_017371",
        "NC_000921",
        "NC_017361",
        "NZ_CP011330",
        "NC_017742",
        "NC_017374",
        "NZ_CP011486",
        "NZ_CP011484",
        "NC_017381",
    ],
    "prochlorococcus_marinus": [],
    "mycobacterium_tuberculosis": [],
    "escherichia_coli": [],
}

# function to extract accession number from names of genbank files in the folder
def get_acc_nums(fld):
    """Given the folder containing the genbank file, extracts the basename"""
    files = glob.glob(fld + "/*.gbk")
    p = re.compile("/?([^/]+)\.gbk$")
    return [re.search(p, f).groups()[0] for f in files]


if __name__ == "__main__":

    # ~~~~~ part 1): select strain pairs for pairwise graphs vs projection comparison  ~~~~~

    # seed random number generator for reproducibility
    rng = np.random.default_rng(0)

    # total number of strains to consider
    # and total number of pairs
    N_tot = 50
    N_pairs = 50

    # collect all accession numbers
    acc_nums = {s: get_acc_nums(f"panx_data/{s}/input_GenBank") for s in species}

    # exclude selected accession number
    for s in species:
        acc_nums[s] = [x for x in acc_nums[s] if x not in exclude[s]]

    # pick random set
    selected_acc_nums = {}
    for s in species:
        acc = acc_nums[s].copy()
        Nsel = min(N_tot, len(acc))
        # using shuffle so that if N_pairs is increased the same strains and pairs
        # are retained. Useful later to avoid re-computing many projections
        rng.shuffle(acc)
        selected_acc_nums[s] = list(sorted(acc[:Nsel]))

    # pick random pairs
    pairs = {}
    for s in species:
        acc = selected_acc_nums[s]
        Na = len(acc)
        N_pick = min(N_pairs, int(Na * (Na - 1) / 2))
        all_pairs = list(itt.combinations(acc, 2))
        rng.shuffle(all_pairs)
        pairs[s] = all_pairs[:N_pick]

    # save to json
    out_fname = "config/projection_strains.json"
    jdict = {s: {"strains": selected_acc_nums[s], "pairs": pairs[s]} for s in species}
    with open(out_fname, "w") as f:
        json.dump(jdict, f, indent=4)

    # ~~~~~ part 2): select strain sets for incremental size pangraph properties  ~~~~~

    species = "escherichia_coli"
    strains = acc_nums[species]  # list of strains
    sizes = [2, 4, 8, 16, 32, 64, 128]  # incremental list of sizes
    n_trials = 10  # n. trials for each size

    # build nested dictionary {size -> n.trial -> [list of strains]}
    jdict = defaultdict(dict)
    for s, t in itt.product(sizes, range(n_trials)):
        # seed random number generator for reproducibility
        # will always have the same results for pairs of size/trial
        seed = hash((s, t))
        if seed < 0:
            seed = -seed
        rng = np.random.default_rng(seed)
        # pick s strains
        jdict[s][t] = list(rng.choice(strains, size=s, replace=False))
    jdict = dict(jdict)

    # save json file
    out_fname = "config/incremental_size_strains.json"
    with open(out_fname, "w") as f:
        json.dump(jdict, f, indent=4)
