from .utils import Pangraph
from .fasta_utils import FastaRecord, write_to_file


def block_to_consensus_records(pangraph: Pangraph) -> list[FastaRecord]:
    records = []
    for n, (block_id, block) in enumerate(pangraph.blocks.items()):
        record_name = f"{block_id}"
        records.append(FastaRecord(name=record_name, seq=block.consensus, idx=n))
    return records


def export_block_consensus(pangraph: Pangraph, file_path: str):
    records = block_to_consensus_records(pangraph)
    write_to_file(records, file_path)
