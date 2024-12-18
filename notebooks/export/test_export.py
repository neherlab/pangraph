import pytest
import json
import tempfile

from .utils import Pangraph
from .export_block_consensus import export_block_consensus
from .export_block_sequences import export_block_seqs
from .export_core_aln import core_block_aln
from .fasta_utils import read_from_file


@pytest.fixture
def load_graph():
    fname = "notebooks/export/test_graph.json"
    with open(fname, "r") as f:
        data = json.load(f)
    graph = Pangraph.from_json_dict(data)
    return graph


def test_export_block_consensus(load_graph):
    graph = load_graph
    with tempfile.TemporaryDirectory() as tmpdirname:
        file_path = f"{tmpdirname}/block_consensus.fa"
        export_block_consensus(graph, file_path)
        records = read_from_file(file_path)
        assert len(records) == len(graph.blocks)
        assert set([int(record.name) for record in records]) == set(graph.blocks.keys())
        for record in records:
            bid = int(record.name)
            assert bid in graph.blocks
            assert record.seq == graph.blocks[bid].consensus


@pytest.mark.parametrize("aligned", [True, False])
def test_export_block_sequences(load_graph, aligned):
    graph = load_graph
    with tempfile.TemporaryDirectory() as tmpdirname:
        export_block_seqs(graph, tmpdirname, aligned)
        for block_id, block in graph.blocks.items():
            records = read_from_file(f"{tmpdirname}/block_{block_id}.fa")
            assert len(records) == len(block.alignments)
            for record in records:
                node_id = int(record.name.split()[0])
                assert node_id in block.alignments

                if aligned:
                    # alignment must have the length of the consensus sequence
                    assert len(record.seq) == len(block.consensus)
                else:
                    # unaligned sequences must have the length of the consensus,
                    # minus the length of deletions, plus the length of insertions
                    seq_len = len(block.consensus)
                    edits = block.alignments[node_id]
                    for D in edits.dels:
                        seq_len -= D.length
                    for I in edits.ins:
                        seq_len += len(I.ins)
                    assert len(record.seq) == seq_len, f"node_id: {node_id}"


@pytest.mark.parametrize("aligned", [True, False])
def test_core_block_aln(load_graph, aligned):
    graph = load_graph
    guide_strain = "pCAV1344-40"
    records = core_block_aln(graph, guide_strain, aligned=aligned)
    assert len(records) == len(graph.paths)
    assert set([record.name for record in records]) == set(
        [path.name for path in graph.paths.values()]
    )
    if aligned:
        L0 = len(records[0].seq)
        for record in records:
            assert len(record.seq) == L0
