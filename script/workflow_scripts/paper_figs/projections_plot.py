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
    parser.add_argument("--csv", help="input stats", nargs="+", type=str)
    parser.add_argument("--pdf", help="output pdf file", type=str)
    return parser.parse_args()


def extract_species(csv):
    name = re.search(r"benchmark/[^/_]+_([^\.]+)\.full\.csv", csv).group(1)
    name = name.split("_")
    sp_name = name[0][0].upper() + ". " + name[1].capitalize()
    return sp_name


def load_dfs(csvs):
    data = {}
    for csv in csvs:
        species = extract_species(csv)
        df = pd.read_csv(csv, index_col=0, header=[0, 1])
        data[species] = df
    return data


def ax_boxplot(df, species, ax):

    sdf = df.reset_index(drop=True).unstack().reset_index()
    sdf = sdf.drop(columns=["level_1"])
    sdf = sdf.rename(columns={"level_0": "kind", 0: "value"})

    # box-plots
    sns.boxplot(ax=ax, data=sdf, x="value", y="kind", color="C0", width=0.6)
    ax.set_ylabel("")
    ax.set_xlabel("")
    ax.grid(axis="x", alpha=0.4)

    # report mean and standard deviation
    # stats = {"mean": df.mean().to_dict(), "std": df.std().to_dict()}
    # M = sdf["value"].max()
    # cols = df.columns
    # for nc, c in enumerate(cols):
    #     mean, std = stats["mean"][c], stats["std"][c]
    #     ax.text(M * 1.1, nc, f"{mean:.2} $\pm$ {std:.2}")


def add_column_header(ax, header):
    ax.annotate(
        header,
        xy=(0.5, 1),
        xytext=(0, 5),
        xycoords="axes fraction",
        textcoords="offset points",
        ha="center",
        va="baseline",
    )


def projection_plot(data, savename):
    NS = len(data)  # number of species

    fig, axs = plt.subplots(NS, 2, figsize=(8, 1.5 * NS), sharey=True, sharex="col")

    for ns, species in enumerate(data):
        df = data[species]

        ax = axs[ns, 0]
        ax_boxplot(df["fraction of total genome length"], species, ax)
        ax.text(0, 0, f"{species}", weight="bold", size="small")

        ax = axs[ns, 1]
        ax_boxplot(df["average segment size (kbp)"], species, ax)
        ax.set_xscale("log")

    add_column_header(axs[0, 0], "fraction of total genome length")
    add_column_header(axs[0, 1], "average segment size (kbp)")

    axs[0, 1].set_xticks(10.0 ** np.arange(-3, 4))

    sns.despine(fig)
    plt.tight_layout()
    plt.savefig(savename)
    plt.close(fig)


if __name__ == "__main__":

    args = parse_args()

    # load data
    data = load_dfs(args.csv)

    # produce plot
    projection_plot(data, args.pdf)
