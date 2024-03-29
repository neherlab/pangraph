import argparse
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

color_core = "#fb8500"
color_all = "#8d99ae"
color_pangraph = "#219EBC"
color_stats = "#023047"


def parse_args():
    parser = argparse.ArgumentParser(
        description="""plot statistics for strain projections"""
    )
    parser.add_argument("--csv", help="input csv file", type=str)
    parser.add_argument("--pdf_all", help="output figure will all stats", type=str)
    parser.add_argument("--pdf_paper", help="output figure for the paper", type=str)
    return parser.parse_args()


def plot_all(df, savename):

    xvar = "n. genomes"
    variables = [
        "avg. fasta file size (Mb)",
        "avg. genome length (bp)",
        "n. core blocks",
        "n. blocks",
        "pangenome length (bp)",
        "avg. block length (bp)",
        "N50 (bp)",
        "L50",
        "F50",
        "fract. core pangenome",
        "len. core pangenome",
        "pangraph file size (Mb)",
        "avg. genome length (Mbp)",
        "filesize compression",
        "sequence compression",
        "pangenome length (Mbp)",
        "fract. core blocks",
        "fract. core genome",
    ]

    NV = len(variables)
    fig, axs = plt.subplots(1, NV, figsize=(NV * 3, 3), sharex=True)

    for nv, yvar in enumerate(variables):

        ax = axs[nv]
        sns.lineplot(
            ax=ax,
            data=df,
            x=xvar,
            y=yvar,
            ci="sd",
            marker=".",
            markeredgecolor=None,
        )
        ax.grid(alpha=0.2)
        ax.set_ylim(bottom=0)

    # tweak second column ax
    ax.set_xscale("log")

    # despine and save
    sns.despine(fig)
    plt.tight_layout()
    plt.savefig(savename)
    plt.close(fig)


def plot_line(ax, df, xvar, yvar, label, color):
    """Utility function to plot the mean and standard
    deviation of a quantity in the dataframe."""
    sdf = df[[xvar, yvar]].copy()
    means = sdf.groupby(xvar).mean()[yvar]
    x = means.index
    stds = sdf.groupby(xvar).std()[yvar].to_numpy()
    ax.plot(x, means, ".-", color=color, label=label)
    ax.fill_between(x, means - stds, means + stds, alpha=0.2, color=color)


def plot_n_blocks_stats(df, ax):

    xvar = "n. genomes"
    yvars = {
        "n. blocks": (color_pangraph, "tot"),
        "n. core blocks": (color_core, "core"),
        "L50": (color_stats, "L50"),
    }
    for yvar in yvars:
        color, label = yvars[yvar]
        plot_line(ax, df, xvar, yvar, label, color)
    ax.grid(alpha=0.2)
    ax.legend()
    ax.set_xlabel("n. genomes")
    ax.set_ylabel("n. pancontigs")
    ax.set_yscale("log")


def plot_block_len_stats(df, ax):

    xvar = "n. genomes"
    yvars = {
        "avg. block length (bp)": (color_pangraph, "avg."),
        "avg. core block length (bp)": (color_core, "avg. core"),
        "N50 (bp)": (color_stats, "N50"),
    }
    for yvar in yvars:
        color, label = yvars[yvar]
        plot_line(ax, df, xvar, yvar, label, color)
    ax.grid(alpha=0.2)
    ax.legend()
    ax.set_xlabel("n. genomes")
    ax.set_ylabel("pancontig size (bp)")
    ax.set_yscale("log")


def plot_pangenome_size(df, ax):

    sdf = df.copy()
    sdf["pang_len"] = sdf["pangenome length (Mbp)"] * 1e6
    sdf["tot_gen_len"] = sdf["avg. genome length (Mbp)"] * 1e6 * sdf["n. genomes"]

    xvar = "n. genomes"
    yvars = {
        "tot_gen_len": (color_all, "all genomes"),
        "pang_len": (color_pangraph, "pangenome"),
        "len. core pangenome": (color_core, "core pangenome"),
    }
    for yvar in yvars:
        color, label = yvars[yvar]
        plot_line(ax, sdf, xvar, yvar, label, color)
    ax.grid(alpha=0.2)
    ax.legend()
    ax.set_xlabel("n. genomes")
    ax.set_ylabel("genome size (bp)")
    ax.set_yscale("log")


def add_panel_label(ax, nax):
    """Add letter to the panel"""
    l = chr(65 + nax)
    ax.text(-0.18, 0.95, l, transform=ax.transAxes, size=15, weight="bold")


def plot_summary(df, savename):

    fig, axs = plt.subplots(1, 3, figsize=(9, 3), sharex=True)

    ax = axs[0]
    plot_n_blocks_stats(df, ax)

    ax = axs[1]
    plot_block_len_stats(df, ax)

    ax = axs[2]
    plot_pangenome_size(df, ax)

    for i in range(3):
        add_panel_label(axs[i], i)

    # tweak second column ax
    ax.set_xscale("log")

    # despine and save
    sns.despine(fig)
    plt.tight_layout()
    plt.savefig(savename)
    plt.close(fig)


if __name__ == "__main__":

    args = parse_args()

    # load dataframe
    df = pd.read_csv(args.csv)

    # produce plot with all statistics
    plot_all(df, args.pdf_all)

    # produce plot for the paper
    plot_summary(df, args.pdf_paper)
