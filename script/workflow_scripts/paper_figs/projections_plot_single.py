import argparse
import re
import json
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


def parse_args():
    parser = argparse.ArgumentParser(
        description="""plot statistics for strain projections"""
    )
    parser.add_argument("--csv", help="input comparison stats", type=str)
    parser.add_argument("--kmer_dist", help="input kmer distance dataframe", type=str)
    parser.add_argument(
        "--core_div",
        help="json with core genome divergence",
        type=str,
    )
    parser.add_argument("--klen", help="kmer size", type=int)
    parser.add_argument("--pdf", help="output pdf file", type=str)
    return parser.parse_args()


def extract_species(csv):
    """Given the csv file name extracts and returns the formatted species name."""
    name = re.search(r"benchmark/[^/_]+_([^\.]+)\.full\.csv", csv).group(1)
    name = name.split("_")
    sp_name = name[0][0].upper() + ". " + name[1].capitalize()
    return sp_name


def load_comparison_stats(csv):
    """Given a filename loads the comparison dataframe and species name"""
    species = extract_species(csv)
    df = pd.read_csv(csv, index_col=0, header=[0, 1])
    return species, df


def add_kmer_distance_to_df(comp_df, kmer_df):
    """Add to the dataframe the (corrected) shared kmer fraction."""
    kmer_dist = kmer_df["corrected_f"].to_dict()
    # kmer_dist = kmer_df["f"].to_dict()
    idx = comp_df.index.to_list()
    sh_fract = []
    for i in idx:
        # extract pair of strains
        s1, s2 = re.search(r"__([^-]+)-([^-]+)\.json\|", i).groups()
        # find corresponding shared kmer fraction
        d = kmer_dist[(s1, f"{s1}-{s2}")]
        sh_fract.append(d)
    comp_df[
        ("fraction of total genome length", "shared kmer fraction (corrected)")
        # ("fraction of total genome length", "shared kmer fraction")
    ] = sh_fract
    return comp_df


def ax_boxplot(df, ax):
    """Boxplot the relevant variables given the sub-dataframe and the ax."""

    # unstack and reformat dataframe
    sdf = df.reset_index(drop=True).unstack().reset_index()
    sdf = sdf.drop(columns=["level_1"])
    sdf = sdf.rename(columns={"level_0": "kind", 0: "value"})

    # perform box-plot
    sns.boxplot(
        ax=ax,
        data=sdf,
        x="value",
        y="kind",
        color="C0",
        width=0.6,
        flierprops={"marker": "."},
    )
    ax.set_ylabel("")
    ax.set_xlabel("")
    ax.grid(axis="x", alpha=0.4)


def projection_plot(comp_df, species, savename):

    # setup figure
    fig, axs = plt.subplots(1, 2, figsize=(8, 2), sharey=False)

    # shared genome fraction
    ax = axs[0]
    sdf = comp_df["fraction of total genome length"]
    ax_boxplot(sdf, ax)

    # add species name
    ax.text(0, 0, f"{species}", weight="bold", size="large")

    # shared segment sizes
    ax = axs[1]
    sdf = comp_df["average segment size (bp)"]
    ax_boxplot(sdf, ax)

    # add titles
    axs[0].set_xlabel("fraction of total genome length")
    axs[1].set_xlabel("average segment size (bp)")

    # tweak axes
    ax.set_xscale("log")

    # despine and save
    sns.despine(fig)
    plt.tight_layout()
    plt.savefig(savename)
    plt.close(fig)


if __name__ == "__main__":

    args = parse_args()

    # load data
    species, comp_df = load_comparison_stats(args.csv)

    # turn kbp into bp
    for k in comp_df.columns:
        if k[0] == "average segment size (kbp)":
            comp_df[("average segment size (bp)", k[1])] = comp_df[k] * 1000

    # load kmer distance
    kmer_df = pd.read_csv(args.kmer_dist, index_col=[0, 1])

    # correct kmer distance for mutations
    with open(args.core_div, "r") as f:
        divergence = json.load(f)
    d = divergence["avg. pairwise divergence"]
    p_correct = np.power(1.0 - d, args.klen)
    kmer_df["corrected_f"] = np.minimum(kmer_df["f"] / p_correct, 1)

    # add shared fraction estimated using kmers
    comp_df = add_kmer_distance_to_df(comp_df, kmer_df)

    # produce plot
    projection_plot(comp_df, species, args.pdf)
