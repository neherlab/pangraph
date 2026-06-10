import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt

from pypangraph.junctions import BackboneJunctions
from pypangraph.plots import linear_junction_plot


def test_linear_junction_plot_smoke(junction_pangraph):
    """Verify linear_junction_plot runs end-to-end on the junction fixture."""
    bj = BackboneJunctions(junction_pangraph, L_thr=500)
    edge = next(iter(bj.edges()))
    fig, ax = plt.subplots()
    linear_junction_plot(ax, bj, edge)
    plt.close(fig)
