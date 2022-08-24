import argparse
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


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
        "n. blocks": ("C0", "all"),
        "n. core blocks": ("C1", "core"),
        "L50": ("C2", "L50"),
    }
    for yvar in yvars:
        color, label = yvars[yvar]
        plot_line(ax, df, xvar, yvar, label, color)
    ax.grid(alpha=0.2)
    ax.legend()
    ax.set_xlabel("n. genomes")
    ax.set_ylabel("n. blocks")
    ax.set_yscale("log")


def plot_block_len_stats(df, ax):

    xvar = "n. genomes"
    yvars = {
        "avg. block length (bp)": ("C0", "avg. block length"),
        "N50 (bp)": ("C1", "N50"),
    }
    for yvar in yvars:
        color, label = yvars[yvar]
        plot_line(ax, df, xvar, yvar, label, color)
    ax.grid(alpha=0.2)
    ax.legend()
    ax.set_xlabel("n. genomes")
    ax.set_ylabel("length (bp)")
    ax.set_yscale("log")


def plot_pangenome_size(df, ax):

    sdf = df.copy()
    sdf["pang_len"] = sdf["pangenome length (Mbp)"] * 1e6
    sdf["tot_gen_len"] = sdf["avg. genome length (Mbp)"] * 1e6 * sdf["n. genomes"]

    xvar = "n. genomes"
    yvars = {
        "tot_gen_len": ("C0", "all genomes"),
        "pang_len": ("C1", "pangenome"),
        "len. core pangenome": ("C2", "core pangenome"),
    }
    for yvar in yvars:
        color, label = yvars[yvar]
        plot_line(ax, sdf, xvar, yvar, label, color)
    ax.grid(alpha=0.2)
    ax.legend()
    ax.set_xlabel("n. genomes")
    ax.set_ylabel("genome size (bp)")
    ax.set_yscale("log")


def plot_summary(df, savename):

    fig, axs = plt.subplots(1, 3, figsize=(9, 3), sharex=True)

    ax = axs[0]
    plot_n_blocks_stats(df, ax)

    ax = axs[1]
    plot_block_len_stats(df, ax)

    ax = axs[2]
    plot_pangenome_size(df, ax)

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
