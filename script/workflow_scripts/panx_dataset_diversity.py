import json
import argparse
import pathlib
import gzip
import numpy as np
import pandas as pd
from Bio import AlignIO
from multiprocessing import Pool
from functools import partial
from itertools import combinations

# Creates a dataframe with the following entries for each gene cluster
# 	"msa" : "GC_0001230_02",

# 	"core" : False,
# 	"duplicated" : False,
# 	"count" : 23,

# 	"len_gap": 1230,
# 	"len_avg": 1143,
# 	"len_nogap": 879,

# 	"avg_kmer_shared_fract": 0.33,
# 	"avg_kmer_jaccardi": 0.29,

# 	"pairwise_divergence": 0.08,
# 	"non_cons_fraction": 0.05,


def parse_args():
    parser = argparse.ArgumentParser(
        description="""Parses the gene cluster file and alignments, and evaluates
        statistics on the diversity of the gene clusters."""
    )
    parser.add_argument("--gc_json", help="input geneCluster.json file", type=str)
    parser.add_argument("--ntot", help="total number of isolates", type=int)
    parser.add_argument("--kmer_l", help="kmer length", type=int)
    parser.add_argument("--gcfolder", help="folder containing gene clusters", type=str)
    parser.add_argument("--n_cores", help="number of cores used", type=int)
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


def parse_geneCluster(gc_file, Ntot):
    """Given the geneCluster.json file and the total number of isolates parses the
    file and returns a list of tuples (gene_cluster_tag, is_core, is_duplicated)"""

    # load gene-cluster file
    with open(gc_file, "r") as f:
        GC = json.load(f)

    # check that Ntot is compatible with the gene cluster counts
    assert_ntot_correct(GC, Ntot)

    # return list of tuples (msa, core, duplicated)
    return [(gc["msa"], gc_is_core(gc, Ntot), gc["dupli"] == "yes") for gc in GC]


def non_cons_fraction(col):
    """Given a column of the alignment returns the fraction of non-consensus sites"""
    _, ct = np.unique(col, return_counts=True)
    n_cons, n_tot = np.max(ct), np.sum(ct)
    f = (n_tot - n_cons) / n_tot
    return f


def pairwise_div(col):
    """Given a column of the alignment returns the pairwise distance of sequences"""
    _, ct = np.unique(col, return_counts=True)
    n_tot = np.sum(ct)
    n_pairs_tot = 0.5 * n_tot * (n_tot - 1)
    n_pairs_same = np.sum(ct * (ct - 1) * 0.5)
    dist = (n_pairs_tot - n_pairs_same) / n_pairs_tot
    return dist


def evaluate_sequence_divergence(aln_mtx):
    """Given the ungapped alignment matrix it evaluates the average non-consensus
    fraction of sites and the average pairwise divergence"""

    # define function to average over alignment functions
    avg_over_aln_columns = lambda f: np.mean(np.apply_along_axis(f, 0, aln_mtx))

    return {
        "non_cons_fraction": avg_over_aln_columns(non_cons_fraction),
        "pairwise_divergence": avg_over_aln_columns(pairwise_div),
    }


def evaluate_shared_kmers(aln, k_len):
    """Given a biopython alignment object and a k-mer length evaluates the average
    shared kmer fraction and the average jaccardi index for all pairs of sequences"""

    # all kmer sets
    ksets = []
    for s in aln:
        seq = str(s.seq.ungap().upper())
        kmer_set = set([seq[l : l + k_len] for l in range(len(seq) - k_len + 1)])
        ksets.append(kmer_set)

    # evaluate shared fraction of kmers
    Fs, Fj = [], []
    for k1, k2 in combinations(ksets, 2):
        ks = len(k1.intersection(k2))
        kt = len(k1.union(k2))
        Fs.append(ks / len(k1))
        Fs.append(ks / len(k2))
        Fj.append(ks / kt)

    return {
        "avg_kmer_shared_fract": np.mean(Fs),
        "avg_kmer_jaccardi": np.mean(Fj),
    }


def alignment_size(aln_obj, aln_mtx):
    """Given the biopython alignment object and the ungapped alignment matrix
    evaluates the alignment depth, the average sequence length, the total
    alignment length and the ungapped alignment length"""

    # sequence lengths
    Ls = [len(str(s.seq.ungap())) for s in aln_obj]

    return {
        "len_gap": aln_obj.get_alignment_length(),
        "len_avg": np.mean(Ls),
        "len_nogap": aln_mtx.shape[1],
        "count": aln_mtx.shape[0],
    }


def analyze_alignment(aln_file, kmer_l):
    """Given an alignment file it parses it and evaluates alignment properties.
    Returns them in a dictionary"""

    # result dictionary
    res = {}

    # load alignment
    with gzip.open(aln_file, "rt") as f:
        aln_obj = AlignIO.read(f, format="fasta")

    # turn into NxL matrix and remove gaps
    aln_mtx = np.array(aln_obj)
    hasgap = lambda x: np.any(x == "-")
    mask = np.apply_along_axis(hasgap, 0, aln_mtx)
    aln_mtx = aln_mtx[:, ~mask]

    # parse lengths and counts
    res |= alignment_size(aln_obj, aln_mtx)

    # if more than one sequence and long enough
    if (res["len_nogap"] > 2 * kmer_l) and (res["count"] > 1):

        res |= evaluate_sequence_divergence(aln_mtx)
        res |= evaluate_shared_kmers(aln_obj, kmer_l)

    return res


def process_single_cluster(gc, kmer_l, fld):
    """Utility function for parallelization. Takes as input a tuple
    gc=(msa, core, dupl), the kmer length and the folder with the alignments.
    Parses the alignment and returns gene cluster properties in a dictionary"""
    msa, core, dupl = gc
    res = {"msa": msa, "core": core, "dupl": dupl}
    aln_fname = fld / f"{msa}_na_aln.fa.gz"
    res |= analyze_alignment(aln_fname, kmer_l)
    return res


def parse_alignments(gc_fld, cluster_list, kmer_l, n_cores):
    """Parses gene clusters in parallel using `n_cores` threads.
    Returns results in a dataframe"""

    # define function to apply to the single cluster
    fld = pathlib.Path(gc_fld)

    with Pool(n_cores) as p:
        cl_df = p.map(
            partial(process_single_cluster, kmer_l=kmer_l, fld=fld), cluster_list, 30
        )

    return pd.DataFrame(cl_df).set_index("msa")


if __name__ == "__main__":

    # parse arguments
    args = parse_args()

    # parse geneCluster.json
    # Extract a list of tuples (genecluster_id, is_core, is_duplcated)
    cluster_list = parse_geneCluster(args.gc_json, args.ntot)

    # parse alignments and evaluate properties (length, sequence divergence...)
    clusters_df = parse_alignments(
        args.gcfolder, cluster_list, args.kmer_l, args.n_cores
    )

    # save dataframe to csv file
    clusters_df.to_csv(args.out)
