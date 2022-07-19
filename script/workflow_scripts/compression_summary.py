# Given a set of fasta files and the corresponding pangenome graph, this script
# evaluates compression statistics. These includes the average genome size, the
# average fasta file size, the number of files, the average pangenome size,

import json
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


def parse_args():
    """Create argument parser and return parsed arguments"""
    parser = argparse.ArgumentParser(
        description="""script to compound compression statistics in a csv table
        and produce a summary figure"""
    )
    parser.add_argument("--jsons", nargs="*", help="input json files")
    parser.add_argument("--csv", help="output csv table")
    parser.add_argument("--pdf", help="output pdf plot")
    return parser.parse_args()


def json_to_dicts(jfile):
    """Given the name of the json file containing the compression statistics,
    extracts the stats and returns them as a dictionary. This dictionary will
    represent a row in the summary dataframe."""

    with open(jfile, "r") as f:
        obj = json.load(f)

    fa_dict = obj["fasta"]
    pg_dict = obj["pangraph"]

    dicts = []
    for kind, val in pg_dict.items():
        val["kind"] = kind
        val.update(fa_dict)
        dicts.append(val)
    return dicts


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


def short_name(x):
    """Abrreviation for the species name"""
    y = x.split("_")
    y[0] = y[0][0]
    y = [str.capitalize(s) for s in y]
    return ". ".join(y)


def barplot(df, ax, y, x="species", hue="kind", logy=True):
    """Utility function for the single barplot"""
    ax.grid(alpha=0.3, zorder=-1)
    sns.barplot(x=x, y=y, hue=hue, data=df, ax=ax)
    if logy:
        ax.set_yscale("log")


def plot(df):
    """Plot summary compression statistics"""
    fig, axs = plt.subplots(2, 2)

    # modify a copy of the dataframe
    df = df.copy()
    df["species"] = df["species"].apply(short_name)

    # create plot arguments
    Kwargs = [
        {"x": "species", "y": "n. genomes", "hue": None},
        {"x": "species", "y": "avg. genome length (Mbp)", "logy": False, "hue": None},
        {"y": "pangenome length (Mbp)", "logy": False},
        {"y": "filesize compression", "logy": False},
        {"y": "sequence compression"},
        {"y": "n. blocks"},
        {"y": "N50 (bp)"},
        {"y": "L50"},
        {"y": "F50"},
        {"y": "fract. core blocks"},
        {"y": "fract. core genome"},
        {"y": "fract. core pangenome"},
    ]

    # create figure
    N = len(Kwargs)
    fig, axs = plt.subplots(N, 1, figsize=(8, N * 4))

    for n, kwargs in enumerate(Kwargs):
        barplot(df=df, ax=axs[n], **kwargs)

    return fig, axs


if __name__ == "__main__":

    # capture arguments
    args = parse_args()

    # compound jsons in df
    df = sum([json_to_dicts(f) for f in args.jsons], [])
    df = pd.DataFrame(df)

    # order by species (n. isolates) and pangraph kind
    species_order = (
        df.groupby("species")["n. genomes"].mean().sort_values().index.to_numpy()
    )
    df["species"] = pd.Categorical(
        df["species"], categories=species_order, ordered=True
    )
    ord_by = ["species", "kind"]
    df = df.sort_values(ord_by)
    col_order = ord_by + [x for x in df.columns if x not in ord_by]
    df = df.reindex(columns=col_order)

    # evaluate additional stats
    evaluate_additional_stats(df)

    # save dataframe
    df.to_csv(args.csv, index=False)

    # produce summary figure
    fig, axs = plot(df)
    plt.tight_layout()
    plt.savefig(args.pdf)
    plt.close(fig)
