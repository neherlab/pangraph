import pytest
import pypangraph as pp
import numpy as np


@pytest.fixture
def load_graph():
    return pp.Pangraph.load_json("tests/data/graph_circ.json")


def test_block_alignment(load_graph):
    pan = load_graph
    bdf = pan.to_blockstats_df()
    core_blocks = bdf[bdf["core"]]
    bid = core_blocks.index[1]
    aln = pan.blocks[bid].get_alignment()
    A = np.array(aln)
    assert A.shape[0] == 15


def test_core_genome_aln(load_graph):
    pan = load_graph
    aln = pan.core_genome_alignment()
    A = np.array(aln)
    assert A.shape[0] == 15
