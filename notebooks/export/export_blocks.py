from .utils import Block, Pangraph, apply_edits_to_ref
from .fasta_utils import FastaRecord, write_to_file
import pathlib


def create_record_id(graph: Pangraph, node_id: str) -> str:
    """
    Given a pangraph and a node id, returns the record id for the export fasta file.
    This is in the format: 'node_id path_name-block_id [start-end|strand]'
    """
    node = graph.nodes[node_id]
    block_id = node.block_id
    path_name = graph.paths[node.path_id].name
    start, end = node.position
    strand = "+" if node.strandedness else "-"
    return f"{node_id} {path_name}-{block_id} [{start}-{end}|{strand}]"


def block_to_aln(
    graph: Pangraph, block: Block, aligned: bool = False
) -> list[FastaRecord]:
    """
    Given a block, returns a list of FastaRecord objects containing a sequence per node.
    If aligned is True, it returns aligned sequences, with gaps for deletions and no insertions.
    If aligned is False, it returns the full unaligned sequences.
    """
    records = []
    for idx, (node_id, edits) in enumerate(block.alignment.items()):
        record_id = create_record_id(graph, node_id)
        seq = apply_edits_to_ref(edits, block.consensus, aligned=aligned)
        records.append(FastaRecord(name=record_id, seq=seq, idx=idx))
    return records


def export_block_seqs(graph: Pangraph, exp_folder: str, aligned: bool):
    """
    Given a pangraph object, exports the block sequences to a fasta file.
    If aligned is True, it exports aligned sequences, with gaps for deletions and no insertions.
    If aligned is False, it exports the full unaligned sequences.
    """
    # create export folder if it does not exist
    fld = pathlib.Path(exp_folder)
    fld.mkdir(parents=True, exist_ok=True)

    # for each block, export the sequence and save the file
    for block_id, block in graph.blocks.items():
        records = block_to_aln(graph, block, aligned=aligned)
        file_path = f"{exp_folder}/block_{block_id}.fa"
        write_to_file(records, file_path)
