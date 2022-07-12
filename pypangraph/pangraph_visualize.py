# The file collects different functions to visualize properties of the pangraph

import matplotlib.pyplot as plt
import numpy as np
import networkx as nx


def block_size_distribution(pan, savename=None):
    """overview of distriution of block sizes. It displays the distribution and
    cumulative distribution of block sizes.
    """

    # list of block lengths
    bl_Ls = [len(block) for block in pan.blocks]

    # logarithmic binning
    ML, mL = np.max(bl_Ls), np.min(bl_Ls)
    log_bins = np.logspace(np.log10(mL), np.log10(ML), 100)

    # create figure
    fig, ax = plt.subplots(2, 1, figsize=(8, 6))

    for axi in ax:
        axi.set_ylabel("block count")
        axi.set_xlabel("block size (bp)")

    # histogram
    axi = ax[0]
    axi.hist(bl_Ls, bins=log_bins)
    axi.set_xscale("log")

    # cumulative distribution
    axi = ax[1]
    axi.hist(bl_Ls, bins=log_bins, density=True, cumulative=True, histtype="step")
    axi.set_xscale("log")
    axi.grid(alpha=0.4)
    axi.text(0.1, 0.9, f"n. blocks = {len(bl_Ls)}", transform=axi.transAxes)
    axi.text(
        0.1,
        0.8,
        f"avg. size = {np.mean(bl_Ls)/1000:.3} $\pm$ {np.std(bl_Ls)/1000:.3} kbp",
        transform=axi.transAxes,
    )
    axi.set_ylabel("cumulative distr.")

    plt.tight_layout()
    if savename is not None:
        plt.savefig(savename, dpi=150, facecolor="w")
    plt.show()


def plot_graph(
    G,
    node_color=None,
    node_color_label=None,
    node_cmap="rainbow",
    node_size=None,
    edge_width=None,
    edge_color="gray",
    strains_overlay=True,
    only_strains=None,
    strains_wiggle=0.1,
    strain_linewidth=4.0,
    strain_alpha=0.2,
    ax=None,
    savename=None,
    show=True,
):
    """General function to plots the pangraph. Takes as main argument a networkx
    object.

    Args:

    - G (networkx Graph) : the graph produced by the pangraph interface

    # optional arguments
    - node_color (list): if a list of number is passed, the nodes are colored
        accordingly. The list must be in the same order as the G.nodes() list.
    - node_color_label (str): the label to add to the colorbar, if nodes are
        colored
    - node_cmap (str or cmpa): the colormap to use for coloring the nodes.
    - node_size (list of int): if specified then node size will be displayed
        accordingly. It must be in the same order as G.nodes().
    - edge_width (list of int): if speciefied then edge width will be drawn
        accordingly. It must be in the same order of G.edges().
    - edge_color (str or list of int): controls the color of edges. Must be in
        the same order as G.edges().
    - strains_overlay (bool): if true then signle transparent paths are drawn on
        top of the graph, one per strain.
    - only_strains (list of str): if present then only the specified strains are
        plotted as overlay
    - strain_linewidth (float): width of strains lines
    - strain_alpha (float): transparency of strain lines
    - strains_wiggle (float): the average displacement of the strains paths from
        the true overlay path. This is done to avoid exactly overlapping all the
        plots.
    - ax (matplotlib ax object): if specified then the plot is drawn on this axis,
        rather than in a new figure.
    - savename (str): if specified then the plot is saved in the specified folder.
    - show (bool): whether to show or not the plot.
    """

    if ax is None:
        fig, ax = plt.subplots(1, 1, figsize=(13, 10))

    # evaluate layout
    pos = nx.kamada_kawai_layout(G)

    # draw edges
    nx.draw_networkx_edges(G, pos, width=edge_width, edge_color=edge_color, ax=ax)

    # draw nodes
    mapp = nx.draw_networkx_nodes(
        G, pos, node_color=node_color, node_size=node_size, ax=ax, cmap=node_cmap
    )
    # optionally plot colorbar
    if node_color is not None:
        plt.colorbar(mapp, shrink=0.2, label=node_color_label, ax=ax)

    # annotate the number of blocks in the graph
    ax.text(0.1, 0.1, f"n. blocks = {len(G.nodes())}", transform=ax.transAxes)

    # draw colored path overlay for each strain
    if strains_overlay:
        i_col = 0
        for strain, path in G.paths.items():

            if only_strains is not None:
                if not strain in only_strains:
                    continue

            # build the path that a strain is following node by node
            x, y = [], []
            for bl in path:
                npos = pos[bl]
                x.append(npos[0])
                y.append(npos[1])
            x.append(x[0])
            y.append(y[0])

            # add noise to avoid exact overlap
            x = np.array(x) + (np.random.rand() - 0.5) * strains_wiggle
            y = np.array(y) + (np.random.rand() - 0.5) * strains_wiggle

            # plot strains overlay
            plt.plot(
                x,
                y,
                linewidth=strain_linewidth,
                alpha=strain_alpha,
                color=f"C{i_col}",
                label=f"{strain}",
            )
            i_col += 1

        ax.legend()

    # remove spines
    for orient in ax.spines:
        ax.spines[orient].set_visible(False)

    plt.tight_layout()
    if savename is not None:
        plt.savefig(savename, dpi=150, facecolor="w")
    if show:
        plt.show()
