from functools import cache

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


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

    # to_sequences() regenerates a block's full alignment; the flanking core blocks
    # recur in every isolate, so memoize per block id to compute each one only once.
    @cache
    def block_sequences(bid):
        """Return the block's node_id -> sequence map, computing it once per block."""
        return pan.blocks[bid].to_sequences()

    records = []
    for iso, junction in edge_map[edge_str].items():
        oriented = junction.to_canonical()

        seq_parts = []
        for ob in oriented.oriented_blocks():
            node_seq = block_sequences(ob.id)[str(ob.node_id)]
            if not ob.strand:
                node_seq = str(Seq(node_seq).reverse_complement())
            seq_parts.append(node_seq)

        record = SeqRecord(Seq("".join(seq_parts)), id=iso, description=edge_str)
        records.append(record)

    return records
