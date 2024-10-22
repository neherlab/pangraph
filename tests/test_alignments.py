import pytest
import pypangraph as pp


@pytest.fixture
def graph():
    fname = "tests/data/plasmid_graph.json"
    pan = pp.Pangraph.load_json(fname)
    return pan


def test_get_alignment(graph):
    block = graph.blocks[68429315432730903]
    aln = block.to_alignment()
    # all sequences have the same length
    assert len(set(len(seq) for seq in aln.values())) == 1
    # all sequences have the same length as the consensus
    assert len(block.consensus()) == len(list(aln.values())[0])
