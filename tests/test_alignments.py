import pytest
import pypangraph as pp
import numpy as np


@pytest.fixture
def graph():
    fname = "tests/data/plasmid_graph.json"
    pan = pp.Pangraph.load_json(fname)
    return pan


def test_get_sequences(graph):
    bid = graph.blocks.keys()[0]
    block = graph.blocks[bid]
    seqs = block.to_sequences()
    assert len(seqs) == 15


def test_get_alignment(graph):
    bid = graph.blocks.keys()[0]
    block = graph.blocks[bid]
    aln = block.to_alignment()
    # all sequences have the same length
    assert len(set(len(seq) for seq in aln.values())) == 1
    # all sequences have the same length as the consensus
    assert len(block.consensus()) == len(list(aln.values())[0])


def test_core_alignment(graph):
    aln = graph.core_genome_alignment()
    A = np.array(aln)
    assert A.shape[0] == 15
    assert A.shape[1] == 64989
