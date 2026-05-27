import pandas as pd


def junctions_dataframe(edge_map, bdf) -> pd.DataFrame:
    """Build a pivot table of accessory length per isolate and edge.

    Args:
        edge_map: dict mapping edge string ID -> dict[isolate, Junction], as cached
            by BackboneJunctions. Only junctions with a flanking edge are present,
            so terminal junctions on linear paths are naturally excluded.
        bdf: Block stats DataFrame with a 'len' column indexed by block id.

    Returns:
        DataFrame with isolates as rows, edge string ids as columns, and the total
        accessory (center) length as values. NaN where an isolate lacks an edge.
        Columns are unsorted; callers may reorder them (e.g. by frequency).
    """
    rows = []
    for edge_str, iso_junctions in edge_map.items():
        for iso, junction in iso_junctions.items():
            length = sum(bdf.loc[node.id, "len"] for node in junction.center.nodes)
            rows.append({"iso": iso, "edge": edge_str, "len": length})
    jdf = pd.DataFrame(rows)
    if jdf.empty:
        return jdf
    return jdf.pivot_table(index="iso", columns="edge", values="len")
