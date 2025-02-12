from .utils import Pangraph
from .fasta_utils import write_to_file, block_to_aln
import pathlib


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
        records = block_to_aln(graph, block, aligned=aligned, record_naming="node")
        file_path = f"{exp_folder}/block_{block_id}.fa"
        write_to_file(records, file_path)


# this function is used by the commands:
# - export block-sequences (aligned=False)
# - export block-alignments (aligned=True)
