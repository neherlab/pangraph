import pandas as pd


def junction_positions(edge_map, pan) -> pd.DataFrame:
    """Find genomic positions of the flanking core blocks for each junction.

    For every (isolate, edge) the flanking blocks are looked up by their unique
    node id, so the coordinates are unambiguous even if a block recurs in a genome.

    The "left" and "right" labels follow each genome's own path order: ``left`` is
    the flank encountered first when walking the path, ``right`` the one after the
    accessory center. The ``strand`` column records whether that path order matches
    the edge's canonical orientation (True) or is inverted relative to it (False).

    Args:
        edge_map: dict mapping edge string ID -> dict[isolate, Junction], as cached
            by BackboneJunctions.
        pan: A Pangraph object (used for node coordinate lookup).

    Returns:
        A DataFrame with MultiIndex (iso, edge) and columns:
        - left_start, left_end: genomic position of the left flanking block
        - right_start, right_end: genomic position of the right flanking block
        - strand: True if the junction is in the canonical edge orientation
    """
    records = []
    for edge_str, iso_junctions in edge_map.items():
        for iso, junction in iso_junctions.items():
            left_row = pan.nodes.df.loc[str(junction.left.node_id)]
            right_row = pan.nodes.df.loc[str(junction.right.node_id)]
            records.append(
                {
                    "iso": iso,
                    "edge": edge_str,
                    "left_start": left_row["start"],
                    "left_end": left_row["end"],
                    "right_start": right_row["start"],
                    "right_end": right_row["end"],
                    "strand": junction.is_canonical(),
                }
            )

    result = pd.DataFrame(records)
    if result.empty:
        return result
    return result.set_index(["iso", "edge"])
