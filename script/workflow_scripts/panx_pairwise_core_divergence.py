import json
import argparse
import pathlib
import gzip
import re
import numpy as np
import pandas as pd
from Bio import AlignIO
from collections import defaultdict


def parse_args():
    parser = argparse.ArgumentParser(
        description="""Parses the gene cluster file and alignments, and evaluates
        statistics on the diversity of the gene clusters."""
    )
    parser.add_argument("--gc_json", help="input geneCluster.json file", type=str)
    parser.add_argument("--ntot", help="total number of isolates", type=int)
    parser.add_argument("--gcfolder", help="folder containing gene clusters", type=str)
    parser.add_argument("--out", help="output csv file", type=str)
    return parser.parse_args()


def gc_is_core(gc, Ntot):
    """Given a gene-cluster and total number of isolates returns a boolean value
    indicating whether the cluster is core"""
    return (gc["dupli"] == "no") and (gc["count"] == Ntot)


def assert_ntot_correct(GC, Ntot):
    """Check that the maximum number of counts for non-duplicated genes is Ntot"""
    counts = [gc["count"] for gc in GC if gc["dupli"] == "no"]
    assert np.max(counts) == Ntot


def parse_core_geneCluster(gc_file, Ntot):
    """Given the geneCluster.json file and the total number of isolates parses the
    file and returns a list of msa names for core gene clusters"""

    # load gene-cluster file
    with open(gc_file, "r") as f:
        GC = json.load(f)

    # check that Ntot is compatible with the gene cluster counts
    assert_ntot_correct(GC, Ntot)

    # return list of tuples (msa, core, duplicated)
    return [gc["msa"] for gc in GC if gc_is_core(gc, Ntot)]


def analyze_alignment(aln):
    """given a biopython alignment object it returns the ungapped length of the
    alignment and a dictionary {(strain1,strain2) : pairwise_divergence} for all
    possible ordered pairs."""

    # turn into NxL matrix and remove gaps
    aln_mtx = np.array(aln)
    hasgap = lambda x: np.any(x == "-")
    mask = np.apply_along_axis(hasgap, 0, aln_mtx)
    aln_mtx = aln_mtx[:, ~mask]

    # length of ungapped alignment
    ungapped_l = aln_mtx.shape[1]

    # list of strains
    strains = [re.search("^([^-]+)-", s.id).group(1) for s in aln]

    # evaluate pairwise divergence dictionary
    pwd_dict = {}
    for ni, si in enumerate(strains):
        divs = (aln_mtx != aln_mtx[ni]).mean(axis=1)
        for nj, sj in enumerate(strains):
            pwd_dict[(si, sj)] = divs[nj]

    return ungapped_l, pwd_dict


def evaluate_pairwise_divergence(gc_fld, cluster_list):
    """Given a list of core gene clusters and the folder containing the alignments,
    it loads each alingment, evaluates pairwise strain diversity, and returns a
    dictionary {(strain1, strain2) : avg pairwise divergence} for all ordered pairs of
    strains."""

    fld = pathlib.Path(gc_fld)

    avg_pw_div = defaultdict(float)
    Ltot = 0
    pairs_set = None

    for gc in cluster_list:

        # load alignment
        aln_fname = fld / f"{gc}_na_aln.fa.gz"
        with gzip.open(aln_fname, "rt") as f:
            aln = AlignIO.read(f, format="fasta")

        # evaluate average ungapped length and pairwise divergence
        l, pw_div = analyze_alignment(aln)

        # check that all pairs are present
        if pairs_set is not None:
            assert pairs_set == set(pw_div.keys()), "Error: pairs sets are not the same"
        else:
            pairs_set = set(pw_div.keys())

        # add tot the total length and divergences
        Ltot += l
        for p in pairs_set:
            avg_pw_div[p] += pw_div[p] * l

    # perform the average and return a dictionary
    avg_pw_div = {p: d / Ltot for p, d in avg_pw_div.items()}

    return avg_pw_div


if __name__ == "__main__":

    # parse arguments
    args = parse_args()

    # parse geneCluster.json and returns a list of core genes msa names
    core_clusters_list = parse_core_geneCluster(args.gc_json, args.ntot)

    # parse core genes alignments and evaluate pairwise divergence for each core gene cluster
    avg_pw_divergence = evaluate_pairwise_divergence(args.gcfolder, core_clusters_list)

    # saves the results as a csv file
    df = pd.DataFrame(
        {"s1": k[0], "s2": k[1], "div": d} for k, d in avg_pw_divergence.items()
    )
    df.to_csv(args.out, index=False)
