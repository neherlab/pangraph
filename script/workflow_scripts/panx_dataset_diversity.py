import json
import argparse
import pathlib
import gzip
import numpy as np
import pandas as pd
from Bio import AlignIO


def parse_args():
    parser = argparse.ArgumentParser(
        description="""parses the geneCluster file from pangraph and return statistics 
        on the core fraction"""
    )
    parser.add_argument("--gc", help="input geneCluster.json file", type=str)
    parser.add_argument("--species", help="species", type=str)
    parser.add_argument("--ntot", help="total number of isolates", type=int)
    parser.add_argument("--kmer_l", help="kmer length", type=int)
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
    clusters = {}
    for gc in GC:
        clusters[gc["msa"]] = {
            "core": gc_is_core(gc, Ntot),
            "duplicated": gc["dupli"] == "yes",
        }

    return clusters


def non_cons_fraction(col):
    """Given a column of the alignment returns the fraction of non-consensus sites"""
    _, ct = np.unique(col, return_counts=True)
    n_cons, n_tot = np.max(ct), np.sum(ct)
    f = (n_tot - n_cons) / n_tot
    return f


def pairwise_div(col):
    """Returns the pairwise distance on the column alignment"""
    _, ct = np.unique(col, return_counts=True)
    n_tot = np.sum(ct)
    n_pairs_tot = 0.5 * n_tot * (n_tot - 1)
    n_pairs_same = np.sum(ct * (ct - 1) * 0.5)
    dist = (n_pairs_tot - n_pairs_same) / n_pairs_tot
    return dist


def pairwise_shared_kmer_fraction(aln, k_len):
    """Given an alignment and k-mer length evaluates the average shared kmer fraction
    and the average jaccardi index for all pairs of sequences"""

    if len(aln) <= 1:
        return {"avg_shared_frac": None, "avg_jaccardi": None}

    ksets = {}
    for s in aln:
        seq = str(s.seq.ungap().upper())
        if len(seq) < 2 * k_len:
            return {"avg_shared_frac": None, "avg_jaccardi": None}
        ksets[s.id] = set([seq[l : l + k_len] for l in range(len(seq) - k_len)])

    Fs = []
    Fj = []
    for id1, k1 in ksets.items():
        for id2, k2 in ksets.items():
            if id1 == id2:
                continue
            ks = len(k1.intersection(k2))
            kt = len(k1.union(k2))
            Fs.append(ks / len(k1))
            Fj.append(ks / kt)

    return {
        "avg_shared_frac": np.mean(Fs),
        "avg_jaccardi": np.mean(Fj),
    }


def analyze_alignment(aln_file, kmer_l):
    """Given an alignment file it parses it and evaluates alignment properties.
    Returns them in a dictionary."""

    # load alignment
    with gzip.open(aln_file, "rt") as f:
        aln = AlignIO.read(f, format="fasta")
    res = pairwise_shared_kmer_fraction(aln, kmer_l)
    aln = np.array(aln)  # shape NxL
    res |= {"len_gap": aln.shape[1], "len_avg": np.sum(aln != "-") / aln.shape[0]}

    # remove columns with gaps
    hasgap = lambda x: np.any(x == "-")
    mask = np.apply_along_axis(hasgap, 0, aln)
    aln = aln[:, ~mask]

    # evaluate length and sequence divergence measures
    res["count"], res["len_nogap"] = aln.shape
    if (res["len_nogap"] < 2 * kmer_l) or (res["count"] <= 1):
        res["non_cons_fraction"] = None
        res["pairwise_divergence"] = None
    else:
        res["non_cons_fraction"] = np.mean(
            np.apply_along_axis(non_cons_fraction, 0, aln)
        )
        res["pairwise_divergence"] = np.mean(np.apply_along_axis(pairwise_div, 0, aln))

    return res


def parse_alignments(gc_fld, clusters, kmer_l):
    """Parse the gene cluster alignments and returns the average diversity and
    average pairwise divergence (length-weighted)."""

    # initialize folder path and containers
    fld = pathlib.Path(gc_fld)

    for n_gc, gc in enumerate(clusters):
        print(f"processing gene cluster {n_gc} / {len(clusters)} : {gc}")
        f = fld / f"{gc}_na_aln.fa.gz"

        # evaluate diversity, pairwise divergence and gene length
        clusters[gc] |= analyze_alignment(f, kmer_l)

    # return average diversity and pairwise divergence, weighted by length
    return clusters


def average_cluster_properties(clusters):

    df = pd.DataFrame(clusters).T
    cdf = df[df["core"]]

    def w_avg(dflab, wlab):
        return np.sum(cdf[dflab] * cdf[wlab]) / np.sum(cdf[wlab])

    return {
        "avg_core_noncons_fraction": w_avg("non_cons_fraction", "len_nogap"),
        "avg_core_pairwise_divergence": w_avg("pairwise_divergence", "len_nogap"),
        "avg_core_shared_kmer_fraction": w_avg("avg_shared_frac", "len_avg"),
        "avg_core_jaccardi_index": w_avg("avg_jaccardi", "len_avg"),
        "pangenome_len": df["len_avg"].sum(),
        "core_genome_len": cdf["len_avg"].sum(),
        "n. genes": len(df),
        "n. core genes": len(cdf),
    }


if __name__ == "__main__":

    # parse arguments
    args = parse_args()

    # initialize results dictionary
    div = {"info": {"species": args.species}}

    # parse gene-cluster file. Obtain a list of core gene clusters and
    clusters = parse_geneCluster(args.gc, args.ntot)

    # parse alignments and evaluate divergence and pairwise diversity
    div["clusters"] = parse_alignments(args.gcfolder, clusters, args.kmer_l)

    # evaluate average properties
    div["info"] = average_cluster_properties(clusters)

    # write to json file
    with open(args.out, "w") as f:
        json.dump(div, f)
