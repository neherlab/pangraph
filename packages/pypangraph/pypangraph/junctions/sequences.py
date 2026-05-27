from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from ..topology_utils import Edge


def junction_sequences(edge_map, pan, edge_str: str) -> list[SeqRecord]:
    """Extract co-oriented sequences spanning a junction.

    For each isolate carrying the given junction, returns a SeqRecord spanning from
    the start of the left core block to the end of the right one. All sequences are
    co-oriented to the canonical edge direction, so the flanking core blocks line up
    across isolates with the accessory content in between. Individual blocks are
    reverse-complemented as needed based on their strand.

    Args:
        edge_map: dict mapping edge string ID -> dict[isolate, Junction], as cached
            by BackboneJunctions.
        pan: A Pangraph object (used for block sequence lookup).
        edge_str: The canonical edge string ID (e.g. "100_f__200_f").

    Returns:
        A list of SeqRecord objects, one per isolate. The record id is the isolate
        name and the description is the edge string ID. Empty if the edge is absent.
    """
    if edge_str not in edge_map:
        return []

    edge = Edge.from_str_id(edge_str)
    records = []
    for iso, junction in edge_map[edge_str].items():
        oriented = junction if junction.is_canonical(edge) else junction.invert()

        all_nodes = [oriented.left] + oriented.center.nodes + [oriented.right]
        seq_parts = []
        for node in all_nodes:
            block = pan.blocks[node.id]
            node_seq = block.to_sequences()[str(node.node_id)]
            if not node.strand:
                node_seq = str(Seq(node_seq).reverse_complement())
            seq_parts.append(node_seq)

        record = SeqRecord(Seq("".join(seq_parts)), id=iso, description=edge_str)
        records.append(record)

    return records
