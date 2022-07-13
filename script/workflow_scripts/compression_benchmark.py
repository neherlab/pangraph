# Given a set of fasta files and the corresponding pangenome graph, this script
# evaluates compression statistics. These includes the average genome size, the
# average fasta file size, the number of files, the average pangenome size,

import os
import re
import json
import argparse
import numpy as np
import pandas as pd
import pypangraph as pp
from Bio import SeqIO


def parse_args():
    """Create argument parser and return parsed arguments"""
    parser = argparse.ArgumentParser(
        description="""
    Evaluates summary compression statistics, comparing a set of fasta files to
    pangenome graphs generated from them. Outputs a json file."""
    )
    parser.add_argument(
        "--fasta", nargs="*", help="fasta files used to build the pangraph"
    )
    parser.add_argument("--pangraphs", nargs="*", help="list of pangenome graphs")
    parser.add_argument("--out_json", help="output json file")
    return parser.parse_args()


def single_fasta_stats(fa):
    """given a fasta file, evaluate filesize and sequence length."""
    # filesize in Mb
    fsize = os.path.getsize(fa) / 1e6
    # genome length in bp
    seq = SeqIO.read(fa, "fasta")
    return {"avg. fasta file size (Mb)": fsize, "avg. genome length (bp)": len(seq)}


def summary_fasta_stats(fas):
    """Given a list of fasta files, evaluates average filesize, genome length
    and number of genomes."""
    stats = [single_fasta_stats(fa) for fa in fas]
    stats = pd.DataFrame(stats).mean().to_dict()
    stats["n. genomes"] = len(fas)
    return stats


def summary_pangraph_stats(pan_files):
    """Given a list of pangraph files evaluates summary statistics and returns them
    in a nested dictionary, one per file. The key is the pangraph kind extracted from
    the filename."""
    stats = {}
    for pf in pan_files:
        # extract pangraph tag from filename
        kind = re.match(".*/pangraph-(.*).json$", pf).group(1)
        # load pangraph and evaluate stats
        pan = pp.Pangraph.load_json(pf)
        stats[kind] = pangraph_stats(pan)
        stats[kind]["pangraph file size (Mb)"] = os.path.getsize(pf) / 1e6
    return stats


def pangraph_stats(pan: pp.Pangraph):
    """Given a pangraph objects, returns a dictionary containing summary statistics."""

    df = pan.to_blockstats_df().sort_values("len", ascending=False)
    stats = {}

    # N. blocks
    # Avg. block length
    # Tot pangenome length
    N = len(df)
    L = df["len"].sum()
    stats["n. blocks"] = N
    stats["pangenome length (bp)"] = L
    stats["avg. block length (bp)"] = df["len"].mean()

    # L50 (n. of contigs that make up 50% pangenome)
    # N50 (len of contig at 50% pangenome length)
    CS = df["len"].cumsum().to_numpy()
    id50 = np.argwhere(CS > L / 2).min()
    stats["N50 (bp)"] = df.iloc[id50]["len"]
    stats["L50"] = id50

    # n. of core blocks
    # fraction of core pangenome
    stats["n. core blocks"] = df["core"].sum()
    stats["fract. core pangenome"] = np.average(df["core"], weights=df["len"])
    stats["len. core pangenome"] = df[df["core"]]["len"].sum()

    return stats


if __name__ == "__main__":

    # capture arguments
    args = parse_args()

    # extract fasta stats
    stats_fa = summary_fasta_stats(args.fasta)

    # extract pangraph stats
    stats_pan = summary_pangraph_stats(args.pangraphs)

    # merge stats
    stats = {"fasta": stats_fa, "pangraph": stats_pan}

    # save to json
    with open(args.out_json, "w") as f:
        json.dump(stats, f, indent=4)

