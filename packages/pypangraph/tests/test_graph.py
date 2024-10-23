import pytest
import pypangraph as pp


def test_load_graph():
    fname = "tests/data/plasmid_graph.json"
    pp.Pangraph.load_json(fname)


@pytest.fixture
def graph():
    fname = "tests/data/plasmid_graph.json"
    pan = pp.Pangraph.load_json(fname)
    return pan


def test_get_strains(graph):
    strains = graph.strains()
    assert len(strains) == 15


def test_blockstats_df(graph):
    df = graph.to_blockstats_df()
    assert df.shape[0] == 137
    assert df["core"].sum() == 27
    assert df["duplicated"].sum() == 10


def test_blockcount_df(graph):
    df = graph.to_blockcount_df()
    assert df.shape[0] == 137
    assert df.shape[1] == 15
    assert df.sum().sum() == 1042


def test_nodes_to_blocks(graph):
    path = graph.paths["RCS49_p1"]
    nodes = path.nodes
    B, S = graph.nodes.nodes_to_blocks(nodes)
    assert len(B) == len(nodes)
    assert len(S) == len(nodes)

    b, s = graph.nodes.node_to_block(nodes[0])
