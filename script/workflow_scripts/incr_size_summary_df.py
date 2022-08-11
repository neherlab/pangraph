import json
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


def parse_args():
    """Create argument parser and return parsed arguments"""
    parser = argparse.ArgumentParser(
        description="script to compound incremental size statistics in a csv table"
    )
    parser.add_argument("--jsons", nargs="+", help="input json files")
    parser.add_argument("--csv", help="output csv dataframe")
    return parser.parse_args()


def json_to_df(jfiles):
    """Given a list of json files, loads them and returns them as a dataframe, each
    file representing a row."""

    df = []
    for jfile in jfiles:
        with open(jfile, "r") as f:
            df.append(json.load(f))
    return pd.DataFrame(df)


def evaluate_additional_stats(df):
    """Evaluates additional statistics (e.g. compression and fraction of core pangenome)
    and adds it to the dataframe."""

    # abbreviations
    NG = df["n. genomes"]
    GL = df["avg. genome length (bp)"]
    PGL = df["pangenome length (bp)"]

    # compression
    df["avg. genome length (Mbp)"] = GL / 1e6
    df["filesize compression"] = df["pangraph file size (Mb)"] / (
        df["avg. fasta file size (Mb)"] * NG
    )
    df["sequence compression"] = PGL / (GL * NG)

    # core genome
    df["pangenome length (Mbp)"] = PGL / 1e6
    df["fract. core blocks"] = df["n. core blocks"] / df["n. blocks"]
    df["fract. core genome"] = df["len. core pangenome"] / GL
    df["fract. core pangenome"] = df["len. core pangenome"] / PGL

    # fraction of blocks
    df["F50"] = df["L50"] / df["n. blocks"]


if __name__ == "__main__":

    # capture arguments
    args = parse_args()

    # compound jsons in df
    df = json_to_df(args.jsons)

    # order by species (n. genomes) and trial number
    ord_by = ["n. genomes", "trial"]
    df = df.sort_values(ord_by)
    col_order = ord_by + [x for x in df.columns if x not in ord_by]
    df = df.reindex(columns=col_order)

    # evaluate additional stats
    evaluate_additional_stats(df)

    # save dataframe
    df.to_csv(args.csv, index=False)
