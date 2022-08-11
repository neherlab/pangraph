import os
import re
import json
import argparse
import numpy as np
import pandas as pd
import pypangraph as pp
from Bio import SeqIO


def parse_args():
    parser = argparse.ArgumentParser(
        description="Evaluates summary statistics for the pangenome graph"
    )
    parser.add_argument("--pangraph", help="input pangenome graph", type=str)
    parser.add_argument("--fasta", help="input fasta files", type=str, nargs="+")
    parser.add_argument("--json", help="output json file", type=str)
    return parser.parse_args()


def fasta_stats(fasta_files):
    """Given a list of fasta files, evaluates average filesize, genome length
    and number of genomes."""
    stats = []
    for fa in fasta_files:
        # filesize in Mb
        fsize = os.path.getsize(fa) / 1e6
        # genome length in bp
        seq = SeqIO.read(fa, "fasta")
        stats.append(
            {"avg. fasta file size (Mb)": fsize, "avg. genome length (bp)": len(seq)}
        )
    # evaluate averages
    stats = pd.DataFrame(stats).mean().to_dict()
    stats["n. genomes"] = len(fasta_files)
    return stats


def pangraph_stats(pan_file):
    """Given a pangraph objects, returns a dictionary containing summary statistics."""

    # load pangraph
    pan = pp.Pangraph.load_json(pan_file)

    # create block stats dataframe
    df = pan.to_blockstats_df().sort_values("len", ascending=False)

    # initialize container
    stats = {}

    # evaluate general statistics
    N = len(df)
    L = df["len"].sum()
    stats["n. blocks"] = N
    stats["pangenome length (bp)"] = int(L)
    stats["avg. block length (bp)"] = df["len"].mean()
    stats["n. strains"] = len(pan.paths)

    # block size statistics:
    # L50 (n. of contigs that make up 50% pangenome)
    # N50 (len of contig at 50% pangenome length)
    CS = df["len"].cumsum().to_numpy()
    id50 = np.argwhere(CS > L / 2).min()
    stats["N50 (bp)"] = int(df.iloc[id50]["len"])
    stats["L50"] = int(id50)

    # core pangenome statistics
    stats["n. core blocks"] = int(df["core"].sum())
    stats["fract. core pangenome"] = np.average(df["core"], weights=df["len"])
    stats["len. core pangenome"] = int(df[df["core"]]["len"].sum())

    # pangraph file size
    stats["pangraph file size (Mb)"] = os.path.getsize(pan_file) / 1e6

    return stats


def extract_trial_number(pangraph_file):
    """Extract size of dataset and trial number from the input filename."""
    size, trial = re.search(r"/(\d+)/(\d+)/pangraph\.json$", pangraph_file).groups()
    return size, trial


if __name__ == "__main__":

    # capture arguments
    args = parse_args()

    # extract fasta stats
    stats_fa = fasta_stats(args.fasta)

    # extract pangraph stats
    stats_pan = pangraph_stats(args.pangraph)

    # merge stats
    stats = stats_fa | stats_pan

    # save trial number and n. strains
    size, trial = extract_trial_number(args.pangraph)
    n_strains = stats["n. strains"]  # from pangraph
    n_genomes = stats["n. genomes"]  # from fasta files
    assert int(size) == n_strains, f"n. strains={n_strains}, size wildcard={size}"
    assert int(size) == n_genomes, f"n. genomes={n_genomes}, size wildcard={size}"
    stats.pop("n. strains", None)  # redundant
    stats["trial"] = int(trial)

    # save to json
    with open(args.json, "w") as f:
        json.dump(stats, f, indent=4, default=str)
