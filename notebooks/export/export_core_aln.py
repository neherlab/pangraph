from .utils import Pangraph
from .fasta_utils import FastaRecord, block_to_aln, concatenate_records, write_to_file


def core_block_aln(
    graph: Pangraph, guide_strain: str, aligned: bool = False
) -> list[FastaRecord]:
    """
    Given a block, returns a list of FastaRecord objects containing a sequence per node.
    If aligned is True, it returns aligned sequences, with gaps for deletions and no insertions.
    If aligned is False, it returns the full unaligned sequences.
    """
    records = []

    core_block_ids = graph.get_core_block_ids()  # id of all core blocks
    guide_path_id = graph.get_path_id(guide_strain)  # id of guide path
    guide_path = graph.paths[guide_path_id]

    guide_block_ids = [
        graph.nodes[node_id].block_id for node_id in guide_path.nodes
    ]  # ids of all blocks in guide path
    guide_block_strand = [
        graph.nodes[node_id].strand for node_id in guide_path.nodes
    ]  # strandedness of all blocks in guide path

    # extract sequences for all core blocks
    for bid, strand in zip(core_block_ids, guide_block_strand):
        if bid not in guide_block_ids:
            continue  # only consider core blocks
        block = graph.blocks[bid]
        block_records = block_to_aln(
            graph, block, aligned=aligned, record_naming="path"
        )
        # reverse-complement if reverse-complemented on guide strain
        if not strand:
            block_records = [record.reverse_complement() for record in block_records]
        records.append(block_records)

    if len(records) == 0:
        # if no record: return empty records
        path_names = graph.get_path_ids()
        return [
            FastaRecord(name=path, seq="", idx=i) for i, path in enumerate(path_names)
        ]
    else:
        # else concatenate them in a single record set
        return concatenate_records(records)


def export_core_alignment(
    graph: Pangraph, guide_strain: str, output_file: str, aligned: bool = False
):
    """
    Export the core block alignment of a pangraph to a fasta file.
    If aligned is True, it exports aligned sequences, with gaps for deletions and no insertions.
    If aligned is False, it exports the full unaligned sequences of core blocks, to be then
    aligned with external tools.
    """
    records = core_block_aln(graph, guide_strain, aligned=aligned)
    write_to_file(records, output_file)
