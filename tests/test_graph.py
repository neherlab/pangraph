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
