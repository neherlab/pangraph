# script to create a summary plot for the performances of pangraph
# as a function of the input dataset size

import argparse
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# parse arguments
def parse_args():
    parser = argparse.ArgumentParser(
        description="""create a summary plot for the performances of pangraph
            as a function of the input dataset size"""
    )
    parser.add_argument("--csv", help="input csv dataframe", type=str)
    parser.add_argument("--pdf", help="output pdf file", type=str)
    return parser.parse_args()


# plots the behavior of a single variable as a function of input dataset size
def single_plot(df, ax, variable):
    g = sns.lineplot(
        data=df,
        x="n-isolates",
        y=variable,
        hue="length",
        ci="sd",
        facet_kws={"legend_out": True},
    )
    ax.set_xscale("log")
    ax.set_xlabel("n. isolates")
    g._legend.set_title("avg. genome length")


# plots the behavior of different variables as a function of input dataset size
def summary_plot(df, savename):

    # variables of interest
    variables = ["mem", "wall-time", "cpu-percent"]
    ylabels = ["max. memory (Gb)", "wall-time (sec)", "cpu-percent"]
    yscales = ["log", "log", None]

    # initialize figure
    fig, axs = plt.subplots(1, 3, figsize=(4 * len(variables), 4))

    # summary plot for each variable
    for nv, var in enumerate(variables):
        ax = axs[nv]
        single_plot(df, ax, var)
        ax.set_yscale(yscales[nv])
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
