import pandas as pd

from .. import topology_utils as tu
from .junction import path_junction_split


def junctions_dataframe(pan, L_thr: int = 500) -> tuple[pd.DataFrame, pd.Series]:
    """Build a junction length DataFrame and edge frequency Series from a pangraph.

    Identifies core backbone blocks (core blocks with length >= L_thr), splits each
    path at these blocks into junctions, and computes the total accessory length in
    each junction.

    Args:
        pan: A Pangraph object.
        L_thr: Minimum block length to be considered part of the backbone (default 500).

    Returns:
        A tuple (jdf, edge_freq) where:
        - jdf: DataFrame with isolates as rows and edges as columns, values are
          the total length of accessory blocks in each junction (NaN if absent).
        - edge_freq: Series with edge string ids as index and occurrence counts
          as values, sorted descending. The columns of jdf are sorted by this order.
    """
    bdf = pan.to_blockstats_df()
    paths = tu.pangraph_to_path_dict(pan)

    def is_core(node_id):
        return (bdf.loc[node_id, "len"] >= L_thr) and bdf.loc[node_id, "core"]

    rows = []
    for iso, path in paths.items():
        junctions = path_junction_split(path, is_core)
        for J in junctions:
            edge = J.flanking_edge()
            if edge is None:
                # terminal junction on a linear path — no flanking edge
                continue
            L = sum(bdf.loc[node.id, "len"] for node in J.center.nodes)
            rows.append({"iso": iso, "edge": edge.to_str_id(), "len": L})
    jdf = pd.DataFrame(rows)
    jdf = jdf.pivot_table(index="iso", columns="edge", values="len")

    # edge frequency: how many isolates have each junction
    edge_freq = jdf.notna().sum(axis=0).sort_values(ascending=False)
    edge_freq.name = "count"

    # sort columns by frequency
    jdf = jdf[edge_freq.index]

    return jdf, edge_freq
