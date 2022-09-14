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
    """Given the csv file name extracts and returns the formatted species name."""
    name = re.search(r"benchmark/[^/_]+_([^\.]+)\.full\.csv", csv).group(1)
    name = name.split("_")
    sp_name = name[0][0].upper() + ". " + name[1].capitalize()
    return sp_name


def load_dfs(csvs):
    """Given the list of filenames loads the dataframes and returns a
    dictionary {species -> df}."""
    data = {}
    for csv in csvs:
        species = extract_species(csv)
        df = pd.read_csv(csv, index_col=0, header=[0, 1])
        data[species] = df
    return data


def ax_boxplot(df, ax):
    """Boxplot the relevant variables given the sub-dataframe and the ax."""

    # unstack and reformat dataframe
    sdf = df.reset_index(drop=True).unstack().reset_index()
    sdf = sdf.drop(columns=["level_1"])
    sdf = sdf.rename(columns={"level_0": "kind", 0: "value"})

    # perform box-plot
    sns.boxplot(ax=ax, data=sdf, x="value", y="kind", color="C0", width=0.6, flierprops={"marker": "."},)
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


# def add_column_title(ax, title):
#     """Add title to column, given the first ax of the column"""
#     ax.annotate(
#         title,
#         xy=(0.5, 1),
#         xytext=(0, 5),
#         xycoords="axes fraction",
#         textcoords="offset points",
#         ha="center",
#         va="baseline",
#     )


def projection_plot(data, savename):

    # setup figure
    NS = len(data)  # number of species
    fig, axs = plt.subplots(NS, 2, figsize=(8, 1.5 * NS), sharey=True, sharex="col")

    # perform the two plots for each species
    for ns, species in enumerate(data):

        df = data[species]  # select species dataframe

        # shared genome fraction
        ax = axs[ns, 0]
        sdf = df["fraction of total genome length"]
        ax_boxplot(sdf, ax)

        # add species name
        ax.text(0, 0, f"{species}", weight="bold", size="small")

        # shared segment sizes
        ax = axs[ns, 1]
        sdf = df["average segment size (kbp)"]
        ax_boxplot(sdf, ax)

    # add titles
    # add_column_title(axs[0, 0], "fraction of total genome length")
    # add_column_title(axs[0, 1], "average segment size (kbp)")
    axs[-1, 0].set_xlabel("fraction of total genome length")
    axs[-1, 1].set_xlabel("average segment size (kbp)")

    # tweak second column ax
    ax.set_xscale("log")
    axs[0, 1].set_xticks(10.0 ** np.arange(-3, 4))

    # despine and save
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
