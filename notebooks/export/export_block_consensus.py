from .utils import Pangraph
from .fasta_utils import FastaRecord, write_to_file


def block_to_consensus_records(pangraph: Pangraph) -> list[FastaRecord]:
    """Given a block it returns a fasta record with the block consensus."""
    records = []
    for n, (block_id, block) in enumerate(pangraph.blocks.items()):
        record_name = f"{block_id}"
        records.append(FastaRecord(name=record_name, seq=block.consensus, idx=n))
    return records


def export_block_consensus(pangraph: Pangraph, file_path: str):
    """
    Given a pangraph object, exports a fasta file with the block consensus sequences.
    """
    records = block_to_consensus_records(pangraph)
    write_to_file(records, file_path)
