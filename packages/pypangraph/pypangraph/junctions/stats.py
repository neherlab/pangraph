from collections import Counter

import pandas as pd

from ..topology_utils import Edge


def _co_oriented_center_paths(edge_str, iso_junctions):
    """Co-orient all center paths for a given edge to canonical direction.

    Junctions for the same edge may appear in opposite orientations on different
    genomes. Before comparing center paths (e.g. for path categories), they must
    be co-oriented so that the same accessory content is recognized regardless of
    genomic strand.

    Args:
        edge_str: Canonical edge string ID.
        iso_junctions: dict mapping isolate name -> Junction.

    Returns:
        dict mapping isolate name -> co-oriented center Walk.
    """
    edge = Edge.from_str_id(edge_str)
    result = {}
    for iso, junction in iso_junctions.items():
        result[iso] = junction.center if junction.is_canonical(edge) else junction.center.invert()
    return result


def _edge_stats(edge_str, iso_junctions, bdf):
    """Compute statistics for a single edge.

    Args:
        edge_str: Canonical edge string ID.
        iso_junctions: dict mapping isolate name -> Junction.
        bdf: Block stats DataFrame (index=block_id, columns include 'len').

    Returns:
        dict with stat column names as keys.
    """
    edge = Edge.from_str_id(edge_str)
    frequency = len(iso_junctions)

    # Co-orient center paths, then group identical ones into categories. The empty
    # path (junction with no accessory blocks) is a category in its own right, so it
    # is counted like any other distinct center path.
    center_paths = _co_oriented_center_paths(edge_str, iso_junctions)
    category_counts = Counter(center_paths.values())

    n_categories = len(category_counts)
    majority_category_freq = max(category_counts.values())
    is_transitive = n_categories == 1
    is_singleton = frequency > 1 and majority_category_freq == frequency - 1

    # Core block lengths
    left_core_length = bdf.loc[edge.left.id, "len"]
    right_core_length = bdf.loc[edge.right.id, "len"]

    # Unique accessory content: collect all distinct block IDs across all isolates
    unique_block_ids = set()
    for junction in iso_junctions.values():
        for node in junction.center.nodes:
            unique_block_ids.add(node.id)
    accessory_length = sum(bdf.loc[bid, "len"] for bid in unique_block_ids)

    return {
        "frequency": frequency,
        "n_categories": n_categories,
        "majority_category_freq": majority_category_freq,
        "is_transitive": is_transitive,
        "is_singleton": is_singleton,
        "left_core_length": left_core_length,
        "right_core_length": right_core_length,
        "accessory_length": accessory_length,
    }


def junction_stats(edge_map, bdf):
    """Compute per-edge statistics for all junctions.

    Args:
        edge_map: dict mapping edge string ID -> dict[isolate, Junction].
            Typically obtained from BackboneJunctions._edge_map.
        bdf: Block stats DataFrame as returned by Pangraph.to_blockstats_df().
            Must have 'len' column and block IDs as index.

    Returns:
        DataFrame with edge string IDs as index and columns:
        - frequency: number of isolates with this junction
        - n_categories: number of distinct center path variants
        - majority_category_freq: count of the most common variant
        - is_transitive: True if only one variant exists
        - is_singleton: True if all but one isolate share the same variant
        - left_core_length: length of left flanking core block
        - right_core_length: length of right flanking core block
        - accessory_length: total unique accessory content (sum of distinct block lengths)

        Sorted by frequency descending.
    """
    # TODO: one could add more stats. E.g. the number of empty paths (no accessory blocks), the average non-empty accessory length.
    records = {}
    for edge_str, iso_junctions in edge_map.items():
        records[edge_str] = _edge_stats(edge_str, iso_junctions, bdf)

    df = pd.DataFrame.from_dict(records, orient="index")
    df.index.name = "edge"
    df = df.sort_values("frequency", ascending=False)

    # Ensure integer types for count columns
    for col in [
        "frequency",
        "n_categories",
        "majority_category_freq",
        "left_core_length",
        "right_core_length",
        "accessory_length",
    ]:
        df[col] = df[col].astype(int)

    return df
