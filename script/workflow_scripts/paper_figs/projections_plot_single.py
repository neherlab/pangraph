import argparse
import re
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
        "--pairwise_div",
        help="csv with core genome pairwise divergence",
        type=str,
    )
    parser.add_argument("--klen", help="kmer size", type=int)
    parser.add_argument("--pdf", help="output pdf file", type=str)
    parser.add_argument("--svg", help="output svg file", type=str)
    return parser.parse_args()


def beautify_species(name):
    """Beautifies the species name"""
    n1, n2 = name.split("_")
    return f"{n1[0].upper()}. {n2}"


def extract_species(csv):
    """Given the csv file name extracts the species name"""
    return re.search(r"benchmark/[^/_]+_([^\.]+)\.full\.csv", csv).group(1)


def load_comparison_stats(csv):
    """Given a filename loads the comparison dataframe and species name"""
    species = extract_species(csv)
    df = pd.read_csv(csv, index_col=0, header=[0, 1])
    return species, df


def add_kmer_distance_to_df(comp_df, kmer_df, pw_div, kmer_l):
    """Add to the dataframe the (corrected) shared kmer fraction."""
    kmer_dist = kmer_df["f"].to_dict()
    # kmer_dist = kmer_df["f"].to_dict()
    idx = comp_df.index.to_list()
    sh_fract = []
    sh_fract_corr = []
    for i in idx:
        # extract pair of strains
        s1, s2 = re.search(r"__([^-]+)-([^-]+)\.json\|", i).groups()

        # find corresponding shared kmer fraction
        d = kmer_dist[(s1, f"{s1}-{s2}")]
        sh_fract.append(d)

        # correct using pairwise divergence
        corr_div = np.power(1 - pw_div[(s1, s2)], kmer_l)
        sh_fract_corr.append(min(d / corr_div, 1))

    main_label = "fraction of total genome length"
    # comp_df[(main_label, "shared kmer fraction")] = sh_fract
    comp_df[(main_label, "shared sequence (kmers)")] = sh_fract_corr
    return comp_df


def ax_boxplot(df, ax):
    """Boxplot the relevant variables given the sub-dataframe and the ax."""

    # unstack and reformat dataframe
    sdf = df.rename(
        columns={
            "agree on sharing": "agree (shared+private)",
            "disagree on sharing": "disagree",
            "shared on both": "agree (shared)",
            "private on both": "agree (private)",
        }
    )
    sdf = sdf.reset_index(drop=True).unstack().reset_index()
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

    # print mean and standard deviation
    means = sdf.groupby("kind")["value"].mean().to_dict()
    stds = sdf.groupby("kind")["value"].std().to_dict()

    K = len(sdf["kind"].unique())
    dK = 1 / K
    for n, c in enumerate(sdf["kind"].unique()):
        txt = f"{means[c]:.3} $\pm$ {stds[c]:.1}"
        txt = scientific_notation(txt)
        ax.text(
            1,
            (1 - dK * 0.6) - n * dK,
            txt,
            # verticalalignment="center",
            transform=ax.transAxes,
        )

    ax.set_ylabel("")
    ax.set_xlabel("")
    ax.grid(axis="x", alpha=0.4)


def scientific_notation(txt):
    def latex_exp(m):
        a1 = "-" if m.group(1) == "-" else ""
        a2 = m.group(2).lstrip("0")
        return r"$\times 10^{" + a1 + a2 + "}$"

    return re.sub("e([+-])([\d][\d])", latex_exp, txt)


def projection_plot(comp_df, species, savename1, savename2):

    # setup figure
    fig, axs = plt.subplots(2, 1, figsize=(7, 5))

    # shared genome fraction
    ax = axs[0]
    sdf = comp_df["fraction of total genome length"]
    ax_boxplot(sdf, ax)

    # add species name
    sp_name = beautify_species(species)
    ax.text(0, 0, f"{sp_name}", weight="bold", size="large")

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
    plt.savefig(savename1)
    plt.savefig(savename2)
    plt.close(fig)


if __name__ == "__main__":

    args = parse_args()

    # load data
    species, comp_df = load_comparison_stats(args.csv)

    # load kmer distance
    kmer_df = pd.read_csv(args.kmer_dist, index_col=[0, 1])

    # add shared fraction estimated using kmers
    pw_div = pd.read_csv(args.pairwise_div, index_col=[0, 1])["div"].to_dict()
    comp_df = add_kmer_distance_to_df(comp_df, kmer_df, pw_div, args.klen)

    # produce plot
    projection_plot(comp_df, species, args.pdf, args.svg)
