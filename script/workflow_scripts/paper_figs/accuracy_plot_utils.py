import json
from math import prod
from tkinter.tix import HList
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from collections import defaultdict
from itertools import product
from matplotlib.lines import Line2D

# load the accuracy dataframe. Optionally keeps only the snps values specified
def load_accuracy_data(json_fname, keep_only_snps=None):
    with open(json_fname, "r") as f:
        data = json.load(f)

    if keep_only_snps is not None:
        data = [x for x in data if x["snps"] in keep_only_snps]

    return data


# creates a dictionary {avg. divergence -> list of costs}
def cost_dictionary(data, mut_factor, keep_only_snps):

    cost_dict = defaultdict(list)
    for d in data:
        snps = d["snps"]
        if not snps in keep_only_snps:
            continue
        costs = d["values"]["costs"]
        divergence = snps * mut_factor
        cost_dict[divergence] += costs

    return dict(cost_dict)


# plots the cumulative distribution of costs, stratified by divergence
def cumulative_cost_plot(costs, ax, legend=True):
    X = sorted(list(costs.keys()))
    Vmax = max([max(v) for v in costs.values()])
    bins = list(np.logspace(-2, np.log10(Vmax) + 0.15, 1000))
    bins = [0] + bins
    cmap = plt.get_cmap("plasma_r")
    N = len(X)

    for n, x in enumerate(X):
        color = cmap(n / (N - 1))
        ax.hist(
            costs[x],
            bins=bins,
            cumulative=True,
            density=True,
            label=f"{x:.3f}",
            histtype="step",
            color=color,
        )
    ax.set_xscale("symlog", linthresh=0.1)
    ax.set_xlim(0, Vmax + 0.1)
    ax.set_ylim(0, 1)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(alpha=0.2)
    if legend:
        ax.legend(
            ncol=2,
            title="avg. divergence",
            fontsize="x-small",
            title_fontsize="small",
            loc="upper left",
        )
    ax.set_xlabel("avg. breakpoint misplacement (bp)")
    ax.set_ylabel("cumul. fraction of strains")


def block_diversity_df(data):
    df = []
    for kernel, k_data in data.items():
        for d in k_data:
            entry = {
                "kernel": kernel,
                "snps": d["snps"],
                "hgt": d["hgt"],
                # "mut_dens": d["values"]["mut_dens"],
                "divergence": d["values"]["dists"],
            }
            df.append(entry)

    return pd.DataFrame(df)


def divergence_vs_snps_rate(df, ax, kernel_title):

    # group dataframe by these values and perform the mean
    gdf = df.groupby(["kernel", "hgt", "snps"]).mean()
    K, H = [df[k].unique() for k in ["kernel", "hgt"]]

    marker = {
        "minimap10": "1",
        "minimap20": "2",
        "mmseqs": "3",
    }

    cmap = plt.get_cmap("plasma")
    color = {h: cmap(nh * 0.95 / (len(H) - 1)) for nh, h in enumerate(H)}

    kwargs = {
        "linewidth": 0.7,
        "ls": ":",
    }

    for k, h in product(K, H):
        # select kernel / hgt pair
        sdf = gdf.loc[k, h].reset_index()
        c, m = color[h], marker[k]
        x = sdf["snps"]
        y = sdf["divergence"]
        ax.plot(x, y, color=c, marker=m, **kwargs)

    # fit line
    sdf = gdf.reset_index()
    mask = sdf["snps"] <= 0.002
    x, y = sdf[mask]["snps"], sdf[mask]["divergence"]
    mut_factor = np.dot(x, y) / np.sum(x**2)
    xmin, xmax = sdf["snps"].min(), sdf["snps"].max()
    xr = np.linspace(xmin, xmax, 10)
    yr = xr * mut_factor
    fit_line = ax.plot(xr, yr, color="gray", ls="--")

    # custom legend
    lines, labels = [], []
    for k in K:
        lines.append(Line2D([0], [0], color="k", marker=marker[k], **kwargs))
        labels.append(kernel_title[k])
    lines.append(fit_line[0])
    labels.append(f"fit (m={mut_factor:.1f})")
    legend1 = ax.legend(
        lines,
        labels,
        fontsize="small",
        title="alignment kernel",
        title_fontsize="small",
        loc="lower right",
    )
    lines, labels = [], []
    for h in H:
        lines.append(Line2D([0], [0], color=color[h], **kwargs))
        labels.append(f"{h:.2f}")
    legend2 = ax.legend(
        lines,
        labels,
        fontsize="small",
        title="hgt rate",
        title_fontsize="small",
        loc="upper left",
    )
    ax.add_artist(legend1)

    # axes setup
    ax.set_ylim(top=sdf["divergence"].max() * 1.1)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(alpha=0.2)
    ax.set_xlabel("simulation mutation rate")
    ax.set_ylabel("avg. pairwise divergence on inferred blocks")

    return mut_factor
