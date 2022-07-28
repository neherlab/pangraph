import json
import argparse
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


def parse_args():
    """Create argument parser and return parsed arguments"""
    parser = argparse.ArgumentParser(
        description="""Collects json files containing the comparison between 
        pairwise graphs and pairwise marginalizations, and produces summary
        statistics."""
    )
    parser.add_argument("--jsons", nargs="+", help="input json files", type=str)
    parser.add_argument("--csv", help="output csv table", type=str)
    parser.add_argument("--pdf", help="output pdf plot", type=str)
    parser.add_argument("--species", help="species for plot title", type=str)
    return parser.parse_args()


def load_stats(jsons):
    """Load the json files and returns a list of data objects"""
    data = []
    for j in jsons:
        with open(j, "r") as f:
            d = json.load(f)
        # forget strain key
        data += list(d.values())
    return data


def data_to_df(data):
    """Turn the nested dictionary data structure into a multi-array dataframe"""
    df = []
    # forget about the strain this was extracted from (remove first key)
    for d in data:
        item = {}
        for k1, v1 in d.items():
            for k2, v in v1.items():
                # key: ( L/N , quantity )
                item[(k2, k1)] = v
        df.append(item)

    # return pd.DataFrame(df)
    return pd.DataFrame(pd.DataFrame(df).to_dict())


def create_summary_df(dfL):
    """Creates a dataframe with fewer columns. These are the ones that will be saved.
    Values represent fraction of genome lengths."""

    tot = dfL["total"]
    sdf = {
        "agree on sharing": dfL["agree on sharing"] / tot,
        "disagree on sharing": dfL["disagree on sharing"] / tot,
        "shared on both": dfL["shared on both"] / tot,
        "private on both": dfL["private on both"] / tot,
        "agree on duplication": dfL["agree on duplication"] / tot,
        "disagree on duplication": dfL["disagree on duplication"] / tot,
        "both are duplicated": dfL["both are duplicated"] / tot,
        "both not duplicated": dfL["both not duplicated"] / tot,
    }
    return pd.DataFrame(sdf)


def plot_summary_df(sdf, savename=None, title=None):
    """Creates two boxplots. One for fraction of shared genome and
    one for fraction of duplicated genome."""

    # prepare sub-df:
    col1 = [
        "agree on sharing",
        "disagree on sharing",
        "shared on both",
        "private on both",
    ]
    col2 = [
        "agree on duplication",
        "disagree on duplication",
        "both are duplicated",
        "both not duplicated",
    ]
    df1 = sdf[col1].unstack().reset_index().drop(columns=["level_1"])
    df2 = sdf[col2].unstack().reset_index().drop(columns=["level_1"])

    mean, std = sdf.mean(), sdf.std()

    fig, axs = plt.subplots(2, 1, sharex=True, figsize=(8, 6))

    # box-plots
    for ax, dfi, col in zip(axs, [df1, df2], [col1, col2]):
        sns.boxplot(ax=ax, data=dfi, x=0, y="level_0", color="C0", width=0.6)
        ax.set_ylabel("")
        ax.set_xlabel("")
        ax.xaxis.grid(True)
        for nc, c in enumerate(col):
            ax.text(1.1, nc, f"{mean[c]:.2} $\pm$ {std[c]:.2}")

    # axes setup
    axs[1].set_xlabel("fraction of genome length")

    st = title.split("_")
    st = f"{st[0].capitalize()[0]}. {st[1].capitalize()}"
    axs[0].set_title(st)

    sns.despine(fig=fig)

    # save and close
    plt.tight_layout()
    plt.savefig(savename)
    plt.close(fig)


if __name__ == "__main__":

    # parse arguments
    args = parse_args()

    # load data and create dataframe
    data = load_stats(args.jsons)
    df = data_to_df(data)

    # create summary dataframe
    dfL = df["L"].copy()
    sdf = create_summary_df(dfL)

    # plot summary of stats
    plot_summary_df(sdf, args.pdf, args.species)

    # save summary
    res = {"mean": sdf.mean(), "std": sdf.std()}
    res = pd.DataFrame(res)
    res.to_csv(args.csv)
