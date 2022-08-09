import json
import argparse
import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict


def parse_args():
    parser = argparse.ArgumentParser(
        description="""create a summary plot for the accuracy of pangraph
            under a different kernel"""
    )
    parser.add_argument("--json", help="input json data", type=str)
    parser.add_argument("--kernel", help="kernel name", type=str)
    parser.add_argument("--pdf", help="output pdf file", type=str)
    parser.add_argument(
        "--mut_factor",
        help="factor to convert snps rate to pairwise diversity",
        type=float,
    )
    parser.add_argument("--snps", help="considered snps values", nargs="+", type=float)
    return parser.parse_args()


# load the accuracy dataframe. Optionally keeps only the snps values specified
def load_accuracy_df(json_fname, keep_only_snps=None):
    with open(json_fname, "r") as f:
        data = json.load(f)

    if keep_only_snps is not None:
        data = [x for x in data if x["snps"] in keep_only_snps]

    return data


def cost_dictionary(data, mut_factor):

    cost_dict = defaultdict(list)
    for d in data:
        snps = d["snps"]
        costs = d["values"]["costs"]
        divergence = snps * mut_factor
        cost_dict[divergence] += costs

    return dict(cost_dict)


def cumulative_cost_plot(costs, ax):
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
    ax.legend(
        ncol=2,
        title="avg. divergence",
        fontsize="x-small",
        title_fontsize="small",
        loc="upper left",
    )
    ax.set_xlabel("avg. breakpoint misplacement (bp)")
    ax.set_ylabel("cumul. fraction of strains")


def accuracy_plot(costs, ker_name, savename):
    fig, ax = plt.subplots(1, 1, figsize=(4.5, 4))

    cumulative_cost_plot(costs, ax)
    title = {
        "minimap10": "minimap2 -asm 10",
        "minimap20": "minimap2 -asm 20",
        "mmseqs": "mmseqs2",
    }
    ax.set_title(title[ker_name])

    plt.tight_layout()
    plt.savefig(savename)
    plt.close(fig)


if __name__ == "__main__":

    # parse arguments
    args = parse_args()

    # load data
    data = load_accuracy_df(args.json, keep_only_snps=args.snps)

    # create cost dictionary
    costs = cost_dictionary(data, args.mut_factor)

    # plot
    accuracy_plot(costs, args.kernel, args.pdf)
