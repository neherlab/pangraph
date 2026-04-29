import pandas as pd

from ..topology_utils import Edge


def _build_node_position_lookup(pan):
    """Build a (block_id, path_id) -> (strand, start, end) lookup dict.

    Core blocks appear exactly once per strain, so the lookup is unambiguous for them.
    """
    lookup = {}
    for _, row in pan.nodes.df.iterrows():
        key = (row["block_id"], row["path_id"])
        lookup[key] = (row["strand"], row["start"], row["end"])
    return lookup


def junction_positions(pan, jdf) -> pd.DataFrame:
    """Find genomic positions of flanking core blocks for each junction.

    For each (isolate, edge) present in jdf, looks up the genomic coordinates of the
    two flanking core blocks and determines the junction strand orientation.

    The "left" and "right" labels refer to the canonical edge orientation: when
    strand=True, the junction appears in the same orientation as the edge definition;
    when strand=False, left/right are swapped relative to the edge definition to
    reflect the actual genomic arrangement.

    The accessory region sits between left_end and right_start, moving in increasing
    coordinate order (wrapping around for circular genomes).

    Args:
        pan: A Pangraph object.
        jdf: Junction length DataFrame as returned by junctions_dataframe (isolates
            as rows, edge string ids as columns).

    Returns:
        A DataFrame with MultiIndex (iso, edge) and columns:
        - left_start, left_end: genomic position of the left flanking block
        - right_start, right_end: genomic position of the right flanking block
        - strand: True if the junction is in the canonical edge orientation
    """
    lookup = _build_node_position_lookup(pan)
    path_ids = {name: path.id for name, path in pan.paths.items()}

    records = []
    for edge_str in jdf.columns:
        edge = Edge.from_str_id(edge_str)
        isolates = jdf[edge_str].dropna().index

        for iso in isolates:
            pid = path_ids[iso]

            # Edge.from_str_id returns string block IDs; the nodes DataFrame
            # stores them as int, so we convert to match.
            left_bid = int(edge.left.id)
            right_bid = int(edge.right.id)

            left_strand, left_start, left_end = lookup[(left_bid, pid)]
            right_strand, right_start, right_end = lookup[(right_bid, pid)]

            # Does the actual genome strand match the edge's node strand?
            same_strand = edge.left.strand == left_strand

            if same_strand:
                records.append({
                    "iso": iso,
                    "edge": edge_str,
                    "left_start": left_start,
                    "left_end": left_end,
                    "right_start": right_start,
                    "right_end": right_end,
                    "strand": True,
                })
            else:
                # Junction is inverted: the edge's "left" block is actually on the
                # right in the genome, and vice versa.
                records.append({
                    "iso": iso,
                    "edge": edge_str,
                    "left_start": right_start,
                    "left_end": right_end,
                    "right_start": left_start,
                    "right_end": left_end,
                    "strand": False,
                })

    result = pd.DataFrame(records)
    if result.empty:
        return result
    return result.set_index(["iso", "edge"])
