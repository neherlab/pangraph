import argparse
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D


def parse_args():
    parser = argparse.ArgumentParser(
        description="""script to produce paper figures for panx benchmark"""
    )
    parser.add_argument("--csv_comp", help="input dataframe for compression")
    parser.add_argument("--csv_summ", help="input dataframe for benchmark")
    parser.add_argument("--pdf_main", help="output pdf plot for main paper")
    parser.add_argument("--pdf_suppl", help="output pdf plot for supplementary")
    return parser.parse_args()


def load_dataframes(csv_comp, csv_summ):
    """Load and concatenate dataframes"""
    dfs = pd.read_csv(csv_summ, index_col=[0, 1])
    dfc = pd.read_csv(csv_comp, index_col=[0, 1])
    df = pd.concat([dfs, dfc], axis=1)
    return df


def short_species_name(name):
    """Return shortened version of species name"""
    sname = name.split("_")
    p1, p2 = sname[0][0].upper(), sname[1].capitalize()
    return f"{p1}. {p2}"


def plot_n_isolates(ax, df, style):
    """Creates bar-plot for the n. of isolates for each species"""
    niso = df["n-isolates"].unstack(level=1).iloc[:, 1].to_dict()
    for i, sp in enumerate(style["sp-order"]):
        ax.bar(i, niso[sp], color="gray")
    set_xlabels(ax, style)
    ax.set_ylabel("n. isolates")


def set_xlabels(ax, style):
    """Sets the x-labels and the y-grid"""
    sord = style["sp-order"]
    ticks = np.arange(len(sord))
    labels = [short_species_name(name) for name in sord]
    ax.set_xticks(ticks)
    ax.set_xticklabels(
        labels,
        rotation=35,
        ha="right",
    )
    ax.set_xlabel("")
    ax.grid(axis="y", alpha=0.2)


def barplot(ax, style, series):
    """Perform single barplot"""

    # column width
    w = 0.15

    # ordered list of species
    S = style["sp-order"]

    for nk, k in enumerate(style["ker-order"]):
        # x-y values
        V = series.loc[:, k].to_dict()
        x = np.arange(len(S)) + (w * nk) - w * 2.5 + 0.05 * (nk > 2)
        y = [V[s] for s in S]
        # color and hatch
        k_type, k_opt = k.split("-")
        color = style["ker-color"][k_type]
        hatch = style["ker-hatch"][k_opt]
        # barplot
        ax.bar(x, y, w, color=color, hatch=hatch, edgecolor="black")

    # setup axis
    set_xlabels(ax, style)


def custom_legend(ax, style, side):
    """Produce double legend: one per alignment kernel and one for options"""

    # capture alignment kernels and options in order
    al_ks, k_opts = [], []
    for k in style["ker-order"]:
        ak, ko = k.split("-")
        if not ak in al_ks:
            al_ks.append(ak)
        if not ko in k_opts:
            k_opts.append(ko)

    # legend for alignment kernels
    elements = []
    for ak in al_ks:
        c = style["ker-color"][ak]
        lab = style["ker-label"][ak]
        elements.append(Patch(facecolor=c, edgecolor=None, label=lab))
    lg1 = ax.legend(
        handles=elements,
        loc=f"upper {side}",
        title="alignment kernel",
        fontsize=8,
        title_fontsize=9,
    )

    # legend for options
    elements = []
    ko_label = {
        "std": r"$\alpha = 100, \beta = 10$",
        "noenergy": r"$\alpha = 0, \beta = 0$",
    }
    for ko in k_opts:
        h = style["ker-hatch"][ko]
        elements.append(
            Patch(facecolor="white", edgecolor="k", label=ko_label[ko], hatch=h)
        )
    lg2 = ax.legend(
        handles=elements,
        loc=f"center {side}",
        title="pseudo-energy",
        fontsize=8,
        title_fontsize=9,
    )
    ax.add_artist(lg1)


def add_panel_label(ax, nax):
    """Add letter to the panel"""
    l = chr(65 + nax)
    ax.text(-0.18, 0.95, l, transform=ax.transAxes, size=15, weight="bold")


def plot_main(df, style, savename):
    """Main text plot"""

    NX, NY = 3, 1
    fig, axs = plt.subplots(
        NY, NX, figsize=(NX * 3.8, NY * 3.3), sharex=True, squeeze=False
    )

    # plot legend
    ax = axs[0, -1]
    custom_legend(ax, style, side="right")
    # plot_n_isolates(ax, df, style)

    # setup for each panel
    yvals = [
        "wall-time",
        "fract. core pangenome",
        "sequence compression",
    ]
    ylabs = {
        "wall-time": "wall time (s)",
    }
    yscale = {
        "wall-time": "log",
        "sequence compression": "log",
    }

    # plot panels
    for i, yv in enumerate(yvals):
        # select axis
        idx = np.unravel_index(i, (NY, NX))
        ax = axs[idx]
        # set y-label
        ylab = ylabs[yv] if yv in ylabs else yv
        ax.set_ylabel(ylab)
        # perform barplot
        series = df[yv]
        barplot(ax, style, series)
        # yscale
        if yv in yscale:
            ax.set_yscale(yscale[yv])

    # setup axes
    sns.despine(fig)
    plt.tight_layout()

    # add panel labels
    for nax in range(NX * NY):
        idx = np.unravel_index(nax, (NY, NX))
        ax = axs[idx]
        add_panel_label(ax, nax)

    plt.savefig(savename)
    plt.close(fig)


def plot_suppl(df, style, savename):
    """Supplementary Information plot"""

    NX, NY = 3, 3
    fig, axs = plt.subplots(NY, NX, figsize=(NX * 4, NY * 3), sharex=True)

    # plot n. isolates in the first panel
    ax = axs[0, 0]
    plot_n_isolates(ax, df, style)
    custom_legend(ax, style, side="left")

    # setup for each panel
    yvals = [
        "wall-time",
        "mem",
        "n. blocks",
        "L50",
        "N50 (bp)",
        "pangenome length (bp)",
        "fract. core pangenome",
        "sequence compression",
    ]
    ylabs = {
        "mem": "max. memory (Gb)",
        "wall-time": "wall time (s)",
    }
    yscale = {
        "wall-time": "log",
        "n. blocks": "log",
        "L50": "log",
        "N50 (bp)": "log",
        "sequence compression": "log",
    }

    # plot panels
    for i, yv in enumerate(yvals):
        # select axis
        idx = np.unravel_index(i + 1, (NY, NX))
        ax = axs[idx]
        # set y-label
        ylab = ylabs[yv] if yv in ylabs else yv
        ax.set_ylabel(ylab)
        # perform barplot
        series = df[yv]
        barplot(ax, style, series)
        # yscale
        if yv in yscale:
            ax.set_yscale(yscale[yv])

    # setup axes
    sns.despine(fig)
    plt.tight_layout()

    # add panel labels
    for nax in range(NX * NY):
        idx = np.unravel_index(nax, (NY, NX))
        ax = axs[idx]
        add_panel_label(ax, nax)

    plt.savefig(savename)
    plt.close(fig)


if __name__ == "__main__":

    # parse arguments
    args = parse_args()

    # load dataframes
    df = load_dataframes(args.csv_comp, args.csv_summ)

    # general style for the plot
    style = {}

    # order of species in the plot (order of n. isolates)
    style["sp-order"] = [
        "prochlorococcus_marinus",
        "mycobacterium_tuberculosis",
        "helicobacter_pylori",
        "klebsiella_pneumoniae",
        "escherichia_coli",
    ]

    # kernel order and style
    style["ker-order"] = [
        "minimap10-std",
        "minimap20-std",
        "mmseqs-std",
        "minimap20-noenergy",
        "mmseqs-noenergy",
    ]
    style["ker-color"] = {
        "minimap10": "#17C3B2",
        "minimap20": "#227C9D",
        "mmseqs": "#A352B7",
    }
    style["ker-hatch"] = {"std": None, "noenergy": "///"}
    style["ker-label"] = {
        "minimap10": "minimap2 asm10",
        "minimap20": "minimap2 asm20",
        "mmseqs": "mmseqs2",
    }

    # plot and save figures
    plot_main(df, style, args.pdf_main)
    plot_suppl(df, style, args.pdf_suppl)

    # df columns
    #     "n-isolates",
    #     "mem",
    #     "wall-time",
    #     "usr-time",
    #     "sys-time",
    #     "cpu-percent",
    #     "command",
    #     "n. blocks",
    #     "pangenome length (bp)",
    #     "avg. block length (bp)",
    #     "N50 (bp)",
    #     "L50",
    #     "n. core blocks",
    #     "fract. core pangenome",
    #     "len. core pangenome",
    #     "pangraph file size (Mb)",
    #     "avg. fasta file size (Mb)",
    #     "avg. genome length (bp)",
    #     "n. genomes",
    #     "avg. genome length (Mbp)",
    #     "filesize compression",
    #     "sequence compression",
    #     "pangenome length (Mbp)",
    #     "fract. core blocks",
    #     "fract. core genome",
    #     "F50",

    # species
    #     'klebsiella_pneumoniae'
    #     'helicobacter_pylori'
    #     'prochlorococcus_marinus'
    #     'mycobacterium_tuberculosis'
    #     'escherichia_coli'

    # kernels
    #      'minimap10-std'
    #      'minimap20-std'
    #      'mmseqs-std'
    #      'minimap20-noenergy'
    #      'mmseqs-noenergy'
