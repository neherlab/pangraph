from collections import defaultdict
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
    """Load the json files and returns a nested dictionary of loaded objects"""
    data = {}
    for j in jsons:
        with open(j, "r") as f:
            d = json.load(f)
        # forget strain key
        k1 = j.split("/")[-1]
        for k2, v in d.items():
            # key: filename.json|strain
            data[f"{k1}|{k2}"] = v
    return data


def data_to_df(data):
    """Turn the nested dictionary data structure into a multi-index dataframe"""
    df, idx = [], []
    # forget about the strain this was extracted from (remove first key)
    for i, d in data.items():
        item = {}
        for k1, v1 in d.items():
            for k2, v in v1.items():
                # key: ( L/N , quantity )
                item[(k2, k1)] = v
        df.append(item)
        idx.append(i)

    # return pd.DataFrame(df)
    return pd.DataFrame(pd.DataFrame(df, index=idx).to_dict())


def create_summary_df(df):
    """Creates a multi-index dataframe with fewer columns. These are the ones that
    will be saved. The first layer of the index corresponds to a different plot"""

    dfL, dfN = df["L"], df["N"]
    k1, k2 = "fraction of total genome length", "average segment size (kbp)"
    vals = [
        "agree on sharing",
        "disagree on sharing",
        "shared on both",
        "private on both",
    ]

    sdf = {}
    for v in vals:
        sdf[(k1, v)] = dfL[v] / dfL["total"]
    for v in vals:
        sdf[(k2, v)] = dfL[v] / (dfN[v] * 1000)

    return pd.DataFrame(sdf)


def plot_summary_df(sdf, savename=None, title=None):
    """Creates a boxplot summary of the dataframe. Each main index corresponds to
    a different plot, and secondary indices are different boxes. For each box the mean
    and standard deviation is reported."""

    # capture mean and std, and keys for main plots
    mean, std = sdf.mean(), sdf.std()
    Ks = sdf.columns.get_level_values(0).unique()

    NK = len(Ks)
    fig, axs = plt.subplots(2, 1, figsize=(8, NK * 3))

    for nk, k in enumerate(Ks):

        # build sub-dataframe
        dfi = sdf[k].unstack().reset_index().drop(columns=["level_1"])
        dfi = dfi.rename(columns={"level_0": "type", 0: "value"})

        # box-plots
        ax = axs[nk]
        sns.boxplot(ax=ax, data=dfi, x="value", y="type", color="C0", width=0.6)
        ax.set_ylabel("")
        ax.set_xlabel(k)
        ax.grid(axis="x", alpha=0.4)

        # report mean and standard deviation
        M = dfi["value"].max()
        cols = sdf[k].columns
        for nc, c in enumerate(cols):
            ax.text(M * 1.1, nc, f"{mean[(k,c)]:.2} $\pm$ {std[(k,c)]:.2}")

    # decorate plot (title and spines)
    st = title.split("_")
    st = f"{st[0].capitalize()[0]}. {st[1].capitalize()}"
    axs[0].set_title(st)
    sns.despine(fig=fig)

    # save and close
    plt.tight_layout()
    plt.savefig(savename)
    plt.close(fig)


def get_top_bottom_N(sdf, N):
    """For each category finds the indices of the top and bottom N values."""
    res = defaultdict(list)
    for n in range(N):
        for k in sdf:
            srt = sdf[k].sort_values()
            srt = srt.dropna()
            res[f"min_{n}"].append(srt.index[n])
            res[f"max_{n}"].append(srt.index[-n - 1])

    idxs = [k for k in sdf]
    res = dict(res)
    for k in res:
        res[k] = pd.Series(data=res[k], index=idxs)
    return res


if __name__ == "__main__":

    # parse arguments
    args = parse_args()

    # load data and create dataframe
    data = load_stats(args.jsons)
    df = data_to_df(data)

    # create summary dataframe
    sdf = create_summary_df(df)

    # plot summary of stats
    plot_summary_df(sdf, args.pdf, args.species)

    # save summary
    res = {"mean": sdf.mean(), "std": sdf.std()} | get_top_bottom_N(sdf, 2)
    res = pd.DataFrame(res)
    res.to_csv(args.csv)
