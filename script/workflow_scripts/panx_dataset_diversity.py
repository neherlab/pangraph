import json
import argparse
import pathlib
import gzip
import numpy as np
from collections import defaultdict
from Bio import AlignIO


def parse_args():
    parser = argparse.ArgumentParser(
        description="""parses the geneCluster file from pangraph and return statistics 
        on the core fraction"""
    )
    parser.add_argument("--gc", help="input geneCluster.json file", type=str)
    parser.add_argument("--species", help="species", type=str)
    parser.add_argument("--ntot", help="total number of isolates", type=int)
    parser.add_argument("--gcfolder", help="folder containing gene clusters", type=str)
    parser.add_argument("--out", help="output json file", type=str)
    return parser.parse_args()


def gc_is_core(gc, Ntot):
    """Given a gene-cluster and total number of entry returns a boolean value
    indicating whether the cluster is core"""
    return (gc["dupli"] == "no") and (gc["count"] == Ntot)


def assert_ntot_correct(GC, Ntot):
    """Check that the maximum number of counts for non-duplicated genes is Ntot"""
    counts = [gc["count"] for gc in GC if gc["dupli"] == "no"]
    assert np.max(counts) == Ntot

def parse_geneCluster(gc_file, Ntot):
    """Given the geneCluster.json file and the total number of isolates
    parses the file and returns the list of core gene-cluster tags and
    a dictionary with the fraction of core genome and pangenome."""

    # load gene-cluster file
    with open(gc_file, "r") as f:
        GC = json.load(f)

    # check that Ntot is compatible with the gene cluster counts
    assert_ntot_correct(GC, Ntot)

    # parse gene clusters
    core_gcs = []
    core_L = defaultdict(int)
    for gc in GC:

        # n. genes
        core_L["n. genes"] += 1

        # length and number of isolates
        l = gc["geneLen"]
        n = gc["count"]

        # total length
        core_L["Ltot"] += l
        core_L["Ltot_n"] += l * n

        if gc_is_core(gc, Ntot):

            # core length
            core_L["Lcore"] += l
            core_L["Lcore_n"] += l * n

            # list of core gene-clusters
            core_gcs.append(gc["msa"])

            # n. genes
            core_L["n. core genes"] += 1

    # evaluate core fractions
    core_F = {} 
    core_F["Fcore_pangenome"] = core_L["Lcore"] / core_L["Ltot"]
    core_F["Fcore_genome"] = core_L["Lcore_n"] / core_L["Ltot_n"]
    core_F["Fcore_genes"] = core_L["n. core genes"] / core_L["n. genes"]

    return core_gcs, core_F


def n_non_consensus(col):
    """Given a column of the alignment evaluates the number of non-consensus sites"""
    _, ct = np.unique(col, return_counts=True)
    n_cons, n_tot = np.max(ct), np.sum(ct)
    n_diff = n_tot - n_cons
    return n_diff


def divergence(aln_file):
    """Given an alignment file evaluates the pairwise divergence on the columns that do
    not contain a gap, and returns it along with the number of ungapped columns."""

    # load alignment
    with gzip.open(aln_file, "rt") as f:
        aln = AlignIO.read(f, format="fasta")
    aln = np.array(aln)  # shape NxL

    # remove gaps
    hasgap = lambda x: np.any(x == "-")
    mask = np.apply_along_axis(hasgap, 0, aln)
    aln = aln[:, ~mask]

    # evaluate length and number of non-consensus sites
    N, l = aln.shape
    N_diff = np.apply_along_axis(n_non_consensus, 0, aln)

    # diversity
    diversity = np.mean(N_diff) / N

    # pairwise divergence
    n_pairs = N * (N - 1) / 2
    pw_divergence = np.mean(N_diff * (N - N_diff) / n_pairs)

    # return diversity, pairwise divergence and length
    return diversity, pw_divergence, l


def parse_alignments(gc_fld, core_gcs):
    """Parse the gene cluster alignments and returns the average diversity and
    average pairwise divergence (length-weighted)."""

    # initialize folder path and containers
    fld = pathlib.Path(gc_fld)
    Diversity, PW_divergence, L = [], [], []

    for n_cg, cg in enumerate(core_gcs):
        print(f"processing core gene cluster {n_cg} / {len(core_gcs)} : {cg}")
        f = fld / f"{cg}_na_aln.fa.gz"

        # evaluate diversity, pairwise divergence and gene length
        d, pwd, l = divergence(f)
        Diversity.append(d)
        PW_divergence.append(pwd)
        L.append(l)

    # return average diversity and pairwise divergence, weighted by length
    return {
        "avg. diversity": np.average(Diversity, weights=L),
        "avg. pairwise divergence": np.average(PW_divergence, weights=L),
    }


if __name__ == "__main__":

    # parse arguments
    args = parse_args()

    print(vars(args))

    # initialize results dictionary
    div = {}
    div["species"] = args.species

    # parse gene-cluster file. Obtain a list of core gene clusters and
    core_gcs, core_fr = parse_geneCluster(args.gc, args.ntot)
    div |= core_fr

    # parse alignments and evaluate divergence and pairwise diversity
    div |= parse_alignments(args.gcfolder, core_gcs)

    # write to json file
    with open(args.out, "w") as f:
        json.dump(div, f)
