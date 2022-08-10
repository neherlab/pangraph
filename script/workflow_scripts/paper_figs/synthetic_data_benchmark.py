# script to create a summary plot for the performances of pangraph
# as a function of the input dataset size

import argparse
import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm


# parse arguments
def parse_args():
    parser = argparse.ArgumentParser(
        description="""create a summary plot for the performances of pangraph
            as a function of the input dataset size"""
    )
    parser.add_argument("--csv", help="input csv dataframe", type=str)
    parser.add_argument("--pdf_main", help="output filename main figure", type=str)
    parser.add_argument(
        "--pdf_suppl", help="output filename supplementary figure", type=str
    )
    return parser.parse_args()


# add a line that goes from the average initial to the average final value of
# the y-axis variable, to display a linear trend.
def plot_lin_trend(df, xlab, ylab, ax):
    xmin = df[xlab].min()
    xmax = df[xlab].max()
    # ymin = df[df[xlab] == xmin][ylab].mean()
    # ymax = df[df[xlab] == xmax][ylab].mean()
    ymax = df[df[xlab] == xmax][ylab].max()

    x = np.logspace(np.log10(xmin), np.log10(xmax), 100)
    # y = ymin + (x - xmin) * (ymax - ymin) / (xmax - xmin)
    y = x * ymax / xmax

    ax.plot(x, y, ls="--", color="gray", label="lin. trend", zorder=-1)


# plots the behavior of a single variable as a function of input dataset size
def single_plot(df, ax, variable, lin_trend=False):

    g = sns.lineplot(
        data=df,
        x="n-isolates",
        y=variable,
        hue="kbp-length",
        hue_norm=LogNorm(),
        ci="sd",
        palette="Blues_d",
        marker=".",
        markeredgecolor=None,
        ax=ax,
    )
    ax.set_xscale("log")
    ax.set_xlabel("n. isolates")
    ax.grid(alpha=0.2)

    # optionally plot a linear trend
    if lin_trend:
        plot_lin_trend(df, "n-isolates", variable, ax)

    g.legend(
        title="avg. isolate length (kbp)",
        fontsize="small",
        title_fontsize="medium",
        ncol=2,
    )


# plots the behavior of different variables as a function of input dataset size
def summary_plot(df, variables, ylabels, yscales, lin_trends, savename):

    # initialize figure
    fig, axs = plt.subplots(
        1, len(variables), figsize=(4 * len(variables), 3.5), squeeze=False
    )

    # summary plot for each variable
    for nv, var in enumerate(variables):

        ax = axs[0, nv]
        single_plot(df, ax, var, lin_trend=lin_trends[nv])

        # customize y-scale and y-label
        ysl = yscales[nv]
        if ysl is not None:
            ax.set_yscale(ysl)
        ax.set_ylabel(ylabels[nv])

    # set axis and
    sns.despine(fig)
    plt.tight_layout()
    plt.savefig(savename)
    plt.close(fig)


if __name__ == "__main__":

    # parse arguments
    args = parse_args()

    # load dataframe
    df = pd.read_csv(args.csv)
    # length in kbp
    df["kbp-length"] = df["length"] // 1000

    # figure version 1
    kwargs = {
        "variables": ["wall-time", "mem", "cpu-percent"],
        "ylabels": ["wall-time (sec)", "max. memory (Gb)", "cpu-percent"],
        "yscales": ["log", "log", None],
        "lin_trends": [True, True, False],
        "savename": args.pdf_suppl,
    }
    summary_plot(df, **kwargs)

    # figure version 2
    kwargs = {
        "variables": ["wall-time"],
        "ylabels": ["wall-time (sec)"],
        "yscales": ["log"],
        "lin_trends": [True],
        "savename": args.pdf_main,
    }
    summary_plot(df, **kwargs)
