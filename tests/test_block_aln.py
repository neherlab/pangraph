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
    block = pan.blocks[bid]
    aln = block.generate_alignment()
    A = np.array(aln)
    assert A.shape[0] == 15


def test_core_genome_aln(load_graph):
    pan = load_graph
    guide_strain = pan.strains()[0]
    aln = pan.core_genome_alignment(guide_strain)
    A = np.array(aln)
    assert A.shape[0] == 15
