import json
from math import prod
from tkinter.tix import HList
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from collections import defaultdict
from itertools import product
from matplotlib.lines import Line2D


def load_accuracy_data(json_fname):
    """Loads the json file containing accuracy data"""
    with open(json_fname, "r") as f:
        data = json.load(f)
    return data


def cost_dictionary(data, conv_factor, keep_only_snps=None, keep_nested_lists=False):
    """Given one of the data dictionaries, returns a dictionary
    { pairwise divergence -> list of avg. breakpoint distances }.
    The pairwise divergence is evaluated by multiplying the snps rate of the
    simulation by a conversion factor. Only values of snp rate that are
    present in `keep_only_snps` are retained.
    If `keep_nested_lists` is set to True then a nested list is returned, each inner
    list corresponding to a different simulation."""

    cost_dict = defaultdict(list)
    # cycle through all simulations
    for d in data:
        # retain only values of snps in the list
        snps = d["snps"]
        if keep_only_snps is not None:
            if not snps in keep_only_snps:
                continue
        # extract list of costs (avg. breakpoint distance per each isolate in the simulation)
        costs = d["values"]["costs"]
        # evaluate divergence from the snps rate
        divergence = snps * conv_factor

        # add the list of costs to the dictionary
        if keep_nested_lists:
            cost_dict[divergence].append(costs)
        else:
            cost_dict[divergence] += costs

    return dict(cost_dict)


def cumulative_cost_plot(costs, ax, legend=True, ylabel=True):
    """Function to plot the cumulative distribution of average breakpoint distance.
    Input:
    - costs: dictionary { divergence -> list of avg. breakpoint distances (~ one per isolate)}
    - ax: ax on which to plot
    - legend: whether to display the legend
    """
    # extract values for the divergence
    D = sorted(list(costs.keys()))
    # maximum value of the cost (~1000)
    cost_max = 1000

    # create the binning
    bins = list(np.logspace(-2, np.log10(cost_max) + 0.15, 1000))
    bins = [0] + bins

    # create dict of colors
    cmap = plt.get_cmap("plasma_r")
    color = {d: cmap(n / (len(D) - 1)) for n, d in enumerate(D)}

    # plot cumulative distributions
    lines, labels = [], []
    for d in D:
        c = color[d]
        ax.hist(
            costs[d],
            bins=bins,
            cumulative=True,
            density=True,
            label=f"{d:.3f}",
            histtype="step",
            color=c,
        )
        lines.append(Line2D([0], [0], color=c))
        labels.append(f"{d:.3f}")
    # plot custom legend
    if legend:
        ax.legend(
            lines,
            labels,
            ncol=1,
            title="avg. divergence",
            fontsize="x-small",
            title_fontsize="small",
            bbox_to_anchor=(1.04, 1),
            loc="upper left",
        )
    # setup axes
    ax.set_xscale("symlog", linthresh=0.1)
    ax.set_xlim(1, cost_max * 1.05)
    ax.set_ylim(0, 1)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(alpha=0.2)
    ax.set_xlabel("avg. breakpoint misplacement (bp)")
    if ylabel:
        ax.set_ylabel("cumul. fraction of isolates")


def block_diversity_df(data):
    """Given the dictionary {alignment kernel -> [list of simulation outcomes]}
    returns a dataframe with columns "kernel", "snps", "hgt", "divergence",
    with one entry per simulation per kernel."""
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


def divergence_vs_snps_rate(df, ax, kernel_title, fit_max_snps):
    """Scatter-plot of the average divergence as a function of the snps rate of the simulation.

    Input:
    - df: block diversity dataframe produced by `block_diversity_df` function.
    - ax: ax on which to plot.
    - kernel_title: dictionary { kernel name -> label in plot }
    - fit_max_snps (float): maximum value of the snps rate to consider for the linear fit.

    Returns:
    - mut_factor: conversion factor from simulation snps rate to average block divergence.
    """

    # group dataframe by kernel, hgt and snps, and evaluate mean divergence
    gdf = df.groupby(["kernel", "hgt", "snps"]).mean()

    # collect all unique values of kernel and hgt rate
    K, H = df["kernel"].unique(), np.sort(df["hgt"].unique())

    # pick marker for each kernel
    marker = {
        "minimap10": "1",
        "minimap20": "2",
        "mmseqs": "3",
    }

    # plot general style
    kwargs = {
        "linewidth": 1.0,
        "ls": ":",
    }

    # dictionary of colors
    cmap = plt.get_cmap("plasma")
    color = {h: cmap(nh * 0.95 / (len(H) - 1)) for nh, h in enumerate(H)}

    # for each pair of kernel and hgt rate, plot the corresponding line
    for k, h in product(K, H):
        # pick dataframe subset
        sdf = gdf.loc[k, h].reset_index()
        # plot the line
        x = sdf["snps"]
        y = sdf["divergence"]
        ax.plot(x, y, color=color[h], marker=marker[k], **kwargs)

    # linear fit of value for small snps rate through the origin
    sdf = gdf.reset_index()
    mask = sdf["snps"] <= fit_max_snps
    x, y = sdf[mask]["snps"], sdf[mask]["divergence"]
    mut_factor = np.dot(x, y) / np.sum(x**2)

    # plot the fit line
    xmin, xmax = sdf["snps"].min(), sdf["snps"].max()
    xr = np.linspace(xmin, xmax, 10)
    yr = xr * mut_factor
    fit_line = ax.plot(xr, yr, color="gray", ls="--", linewidth=1.0)

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
    ax.legend(
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


def mean_and_stderr(values):
    """Given a list of values evaluates the mean and standard error of the mean"""
    mean = np.mean(values)
    stderr = np.std(values) / np.sqrt(len(values))
    return mean, stderr


def median_plot_util(costs, ax, label, color):
    """Given the dictionary of {divergence -> [[avg. isolate displacement] per simulation]}
    plots a line that quantifies the median breakpoint displacement. The median is
    calculated for each simulation and then the mean and standard error are
    reported in the plot."""
    # values of divergence
    divs = list(sorted(costs.keys()))
    mean, stderr = [], []
    # for each value of the divergence
    for d in divs:
        # for each simulation evaluate median displacement
        values = [np.median(c) for c in costs[d] if len(c) > 0]
        # extract mean and stderr over the simulations
        m, s = mean_and_stderr(values)
        mean.append(m)
        stderr.append(s)
    M, S = np.array(mean), np.array(stderr)
    # plot and set y-axis label
    ax.plot(divs, M, ".:", label=label, color=color)
    ax.fill_between(divs, M - S, M + S, alpha=0.2, color=color)
    # medians = [np.median(np.concatenate(costs[d])) for d in divs]
    ax.set_ylabel("median breakpoint misplacement (<1kbp)")


def fraction_plot_util(costs, ax, label, color):
    """Given the dictionary of {divergence -> [[avg. isolate displacement] per simulation]}
    plots a line that quantifies the fraction of costs that are >100bp. This fraction is
    calculated for each simulation and then the mean and standard error are
    reported in the plot."""
    # values of divergence
    divs = list(sorted(costs.keys()))
    mean, stderr = [], []
    # for each value of the divergence
    for d in divs:
        # for each simulation evaluate fraction of average displacements > 100bp
        values = [np.mean(np.array(c) > 100) for c in costs[d] if len(c) > 0]
        # extract mean and stderr over the simulations
        m, s = mean_and_stderr(values)
        mean.append(m)
        stderr.append(s)
    M, S = np.array(mean), np.array(stderr)
    # plot and set y-axis label
    ax.plot(divs, M, ".:", label=label, color=color)
    ax.fill_between(divs, M - S, M + S, alpha=0.2, color=color)
    ax.set_ylabel("fraction of misplaced (>100bp) breakpoints")


def misplacement_vs_divergence(costs, ax, stat):
    """Given the nested dictionary
    {kernel -> divergence -> [[avg. isolate displacement] per simulation]}
    plots either the median displacement (stat='median') or the fraction
    of isolates with average displacement >100bp (stat='fraction').
    """

    legend = {
        "minimap10": "minimap2 asm10",
        "minimap20": "minimap2 asm20",
        "mmseqs": "mmseqs2",
    }
    color = {
        "minimap10": "#17C3B2",
        "minimap20": "#227C9D",
        "mmseqs": "#A352B7",
    }

    # for every kernel
    for k, kernel_costs in costs.items():
        # plot the selected statistic
        if stat == "median":
            median_plot_util(kernel_costs, ax, legend[k], color[k])
        elif stat == "fraction":
            fraction_plot_util(kernel_costs, ax, legend[k], color[k])
        else:
            raise ValueError("stat must be either median or fraction")

    # setup axes
    ax.legend(title="alignment kernel", loc="best")
    # ax.set_xscale("log")
    # ax.set_yscale("log")
    ax.set_xlabel("avg. pairwise divergence")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(alpha=0.2)
