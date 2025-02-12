import pytest
import json
import tempfile

from .utils import Pangraph
from .export_block_consensus import export_block_consensus
from .export_block_sequences import export_block_seqs
from .export_core_aln import core_block_aln, export_core_alignment
from .export_gfa import graph_to_gfa, export_graph_to_gfa
from .fasta_utils import read_from_file


@pytest.fixture
def load_graph():
    fname = "notebooks/export/test_data/test_graph.json"
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


@pytest.mark.parametrize("aligned", [True, False])
def test_export_core_alignment(load_graph, aligned):
    graph = load_graph
    guide_strain = "pCAV1344-40"
    with tempfile.TemporaryDirectory() as tmpdirname:
        output_file = f"{tmpdirname}/core_alignment.fa"
        export_core_alignment(graph, guide_strain, output_file, aligned=aligned)
        records = read_from_file(output_file)
        assert len(records) == len(graph.paths)
        assert set([record.name for record in records]) == set(
            [path.name for path in graph.paths.values()]
        )


@pytest.mark.parametrize("min_len", [None, 1000])
@pytest.mark.parametrize("min_depth", [None, 2])
@pytest.mark.parametrize("export_duplicated", [True, False])
def test_graph_to_gfa(load_graph, min_len, export_duplicated, min_depth):
    gfa = graph_to_gfa(
        load_graph,
        min_len=min_len,
        min_depth=min_depth,
        export_duplicated=export_duplicated,
    )

    # how many blocks after filtering?
    n_blocks = 0
    for block_id, block in load_graph.blocks.items():
        isolates = [
            load_graph.nodes[node_id].path_id for node_id in block.alignments.keys()
        ]
        is_duplicated = len(set(isolates)) < len(isolates)
        if min_len is not None and len(block.consensus) < min_len:
            continue
        if min_depth is not None and len(block.alignments) < min_depth:
            continue
        if not export_duplicated and is_duplicated:
            continue
        n_blocks += 1

    assert len(gfa.segments) == n_blocks


@pytest.mark.parametrize("export_sequence", [True, False])
def test_export_graph_to_gfa(load_graph, export_sequence):
    with tempfile.TemporaryDirectory() as tmpdirname:
        gfa_file = f"{tmpdirname}/graph.gfa"
        export_graph_to_gfa(
            load_graph,
            gfa_file,
            min_len=100,
            min_depth=None,
            export_duplicated=True,
            export_sequence=export_sequence,
        )

        with open(gfa_file, "r") as f:
            lines = f.readlines()

        # check that the number of segments is correct
        n_segments = 0
        for line in lines:
            if line.startswith("S"):
                n_segments += 1
        assert n_segments == len(load_graph.blocks)
