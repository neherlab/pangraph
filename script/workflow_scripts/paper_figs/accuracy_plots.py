import argparse
import matplotlib.pyplot as plt

from accuracy_plot_utils import (
    load_accuracy_data,
    cost_dictionary,
    cumulative_cost_plot,
    block_diversity_df,
    divergence_vs_snps_rate,
    misplacement_vs_divergence,
)


def parse_args():
    parser = argparse.ArgumentParser(
        description="""create a summary plot for the accuracy of pangraph
            under different kernels"""
    )
    parser.add_argument("--mm10", help="input json for minimap10 kernel", type=str)
    parser.add_argument("--mm20", help="input json for minimap20 kernel", type=str)
    parser.add_argument("--mmsq", help="input json for mmseqs2 kernel", type=str)
    parser.add_argument(
        "--snps_suppl",
        help="snps values to consider (supplementary figure)",
        nargs="+",
        type=float,
    )
    parser.add_argument(
        "--snps_main",
        help="snps values to consider (main figure)",
        nargs="+",
        type=float,
    )
    parser.add_argument(
        "--pdf_acc", help="output supplementary accuracy figure", type=str
    )
    parser.add_argument(
        "--pdf_snps",
        help="output supplementary snps vs divergence figure",
        type=str,
    )
    parser.add_argument(
        "--pdf_med",
        help="output median misplacement vs divergence figure",
        type=str,
    )
    parser.add_argument(
        "--pdf_frac",
        help="output misplaced fraction vs divergence figure",
        type=str,
    )
    return parser.parse_args()


def single_accuracy_plot(costs, title, savename):
    """cumulative distribution of breakpoint distance vs sequence divergence
    for a particular alignment kernel."""
    fig, ax = plt.subplots(1, 1, figsize=(4.5, 4))
    cumulative_cost_plot(costs, ax)
    ax.set_title(title)
    plt.tight_layout()
    plt.savefig(savename)
    plt.close(fig)


def comparison_accuracy_plot(costs, titles, savename):
    """comparison of cumulative distributions of breakpoint distance vs sequence divergence
    for all alignment kernel."""
    Nc = len(costs)
    fig, axs = plt.subplots(1, Nc, figsize=(Nc * 3.3, 3), sharey=True)
    for nc, k in enumerate(costs):
        ax = axs[nc]
        cumulative_cost_plot(costs[k], ax, legend=nc == Nc - 1, ylabel=nc == 0)
        ax.set_title(titles[k])
    plt.tight_layout()
    plt.savefig(savename)
    plt.close(fig)


def snps_rate_vs_divergence_plot(df, savename, kernel_title, fit_max_snps):
    """Plot of simulation snps rate vs sequence divergence. Used to find the
    conversion factor between simulation snps rate and sequence divergence."""
    fig, ax = plt.subplots(1, 1, figsize=(5, 4))
    mut_factor = divergence_vs_snps_rate(df, ax, kernel_title, fit_max_snps)
    plt.tight_layout()
    plt.savefig(savename)
    plt.close(fig)
    return mut_factor


def misplacement_plot(costs, savename, stat):
    """Plot of divergence vs:
    - (stat='median') -> median breakpoint displacement
    - (stat='fraction') -> ~ fraction of displaced breakpoints (>100bp)
    """
    fig, ax = plt.subplots(1, 1, figsize=(5, 4))
    misplacement_vs_divergence(costs, ax, stat=stat)
    plt.tight_layout()
    plt.savefig(savename)
    plt.close(fig)


if __name__ == "__main__":

    args = parse_args()

    # extract filenames
    fname = {
        "minimap10": args.mm10,
        "minimap20": args.mm20,
        "mmseqs": args.mmsq,
    }

    titles = {
        "minimap10": "minimap2 asm10",
        "minimap20": "minimap2 asm20",
        "mmseqs": "mmseqs2",
    }

    # load data
    data = {k: load_accuracy_data(v) for k, v in fname.items()}

    # extract block diversity dataframe
    df = block_diversity_df(data)

    # compare snps rate and divergence. Extract conversion factor
    conv_factor = snps_rate_vs_divergence_plot(
        df,
        args.pdf_snps,
        kernel_title=titles,
        fit_max_snps=0.002,
    )

    # extract breakpoint misplacement distance, keep only relevant distance,
    # stratify by sequence divergence
    # creates a dictionary {kernel -> divergence -> [avg. displacement over isolates]}
    costs = {
        k: cost_dictionary(
            v,
            keep_only_snps=args.snps_suppl,
            conv_factor=conv_factor,
            keep_nested_lists=False,
        )
        for k, v in data.items()
    }

    # cumulative distribution of breakpoint misplacement vs sequence divergence
    comparison_accuracy_plot(costs, titles, args.pdf_acc)

    # re-evaluate breakpoint misplacement distance, keep specified snps, leave nested lists
    # for mean/stderr evaluation
    # creates a dictionary {kernel -> divergence -> [[average displacements over isolates] per simulation ]}
    costs = {
        k: cost_dictionary(
            v,
            keep_only_snps=args.snps_main,
            conv_factor=conv_factor,
            keep_nested_lists=True,
        )
        for k, v in data.items()
    }

    # median misplacement vs divergence
    misplacement_plot(costs, args.pdf_med, stat="median")

    # fraction of <100 bp misplacement vs divergence
    misplacement_plot(costs, args.pdf_frac, stat="fraction")
