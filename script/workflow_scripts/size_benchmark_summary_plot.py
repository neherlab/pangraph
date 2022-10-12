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
    parser.add_argument("--pdf", help="output pdf file", type=str)
    return parser.parse_args()

# add a line that goes from the average initial to the average final value of 
# the y-axis variable, to display a linear trend.
def plot_lin_trend(df, xlab, ylab, ax):
    xmin = df[xlab].min()
    xmax = df[xlab].max()
    ymin = df[df[xlab] == xmin][ylab].mean()
    ymax = df[df[xlab] == xmax][ylab].mean()

    x = np.logspace(np.log10(xmin), np.log10(xmax), 100)
    y = ymin + (x - xmin) * (ymax - ymin) / (xmax - xmin)

    ax.plot(x, y, ls="--", color="k", label="lin. trend")


# plots the behavior of a single variable as a function of input dataset size
def single_plot(df, ax, variable, lin_trend=False):
    g = sns.lineplot(
        data=df,
        x="n-isolates",
        y=variable,
        hue="length",
        hue_norm=LogNorm(),
        ci="sd",
        palette="Blues_d",
        marker=".",
        ax=ax,
    )
    ax.set_xscale("log")
    ax.set_xlabel("n. isolates")
    ax.grid(alpha=0.2)

    if lin_trend:
        plot_lin_trend(df, "n-isolates", variable, ax)

    g.legend(title="avg. seq. length", fontsize="small", title_fontsize="medium")


# plots the behavior of different variables as a function of input dataset size
def summary_plot(df, savename):

    # variables of interest
    variables = ["wall-time", "mem", "cpu-percent"]
    ylabels = ["wall-time (sec)", "max. memory (Gb)", "cpu-percent"]
    yscales = ["log", "log", None]

    # initialize figure
    fig, axs = plt.subplots(1, 3, figsize=(4 * len(variables), 4))

    # summary plot for each variable
    for nv, var in enumerate(variables):

        ax = axs[nv]
        lin_trend = var == "wall-time"
        single_plot(df, ax, var, lin_trend=lin_trend)

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

    # plot results and save the figure
    summary_plot(df, args.pdf)
