import pypangraph as pp
import pypangraph.msu as msu
from collections import defaultdict
import pytest


class TestNode:
    def test_eq(self):
        n1 = msu.Node("A", True)
        n2 = msu.Node("A", True)
        n3 = msu.Node("A", False)
        n4 = msu.Node("B", True)
        assert n1 == n2
        assert n1 != n3
        assert n1 != n4

    def test_invert(self):
        n = msu.Node("A", True)
        assert n.invert() == msu.Node("A", False)


class TestPath:
    def test_eq(self):
        n1 = msu.Node("A", True)
        n2 = msu.Node("B", True)
        n3 = msu.Node("C", True)
        p1 = msu.Path([n1, n2, n3], circular=True)
        p2 = msu.Path([n1, n2, n3], circular=True)
        p3 = msu.Path([n1, ~n2, n3], circular=True)
        p4 = msu.Path([n1, n2], circular=True)
        assert p1 == p2
        assert p1 != p3
        assert p1 != p4

    def test_invert(self):
        n1 = msu.Node("A", True)
        n2 = msu.Node("B", True)
        n3 = msu.Node("C", True)
        p = msu.Path([n1, n2, n3], circular=True)
        assert ~p == msu.Path([~n3, ~n2, ~n1], circular=True)

    def test_rotate_to(self):
        n1 = msu.Node("A", True)
        n2 = msu.Node("B", True)
        n3 = msu.Node("C", True)
        n4 = msu.Node("D", True)
        p = msu.Path([n1, n2, n3, n4], circular=True)
        assert p.rotate_to("B", True) == msu.Path([n2, n3, n4, n1], circular=True)
        assert p.rotate_to("D", True) == msu.Path([n4, n1, n2, n3], circular=True)
        assert p.rotate_to("B", False) == msu.Path([~n2, ~n1, ~n4, ~n3], circular=True)
        assert p.rotate_to("D", False) == msu.Path([~n4, ~n3, ~n2, ~n1], circular=True)

    def test_rename_bids(self):
        n1 = msu.Node("A", True)
        n2 = msu.Node("B", False)
        n3 = msu.Node("C", True)
        p = msu.Path([n1, n2, n3], circular=True)
        assert p.rename_bids({"A": "X", "B": "Y", "C": "Z"}) == msu.Path(
            [msu.Node("X", True), msu.Node("Y", False), msu.Node("Z", True)],
            circular=True,
        )


class TestEdge:
    def test_invert(self):
        n1 = msu.Node("A", True)
        n2 = msu.Node("B", True)
        e = msu.Edge(n1, n2)
        assert e.invert() == msu.Edge(n2.invert(), n1.invert())

    def test_eq(self):
        n1 = msu.Node("A", True)
        n2 = msu.Node("B", True)
        e1 = msu.Edge(n1, n2)
        e2 = msu.Edge(n1, n2)
        e3 = msu.Edge(~n1, ~n2)
        e4 = msu.Edge(n1, ~n2)
        e5 = msu.Edge(~n2, ~n1)
        assert e1 == e2
        assert e1 != e3
        assert e1 != e4
        assert e1 == e5


@pytest.fixture
def generate_core_paths():
    A = msu.Node("A", True)
    B = msu.Node("B", True)
    C = msu.Node("C", True)  # invert
    D = msu.Node("D", True)  # invert
    E = msu.Node("E", True)
    F = msu.Node("F", True)
    G = msu.Node("G", True)
    H = msu.Node("H", True)  # invert
    J = msu.Node("J", True)

    p1 = msu.Path([A, B, C, D, E, F, G, H, J], circular=True)
    p2 = msu.Path([A, B, C, D, E, F, G, H, J], circular=True)
    p3 = msu.Path([A, B, ~D, ~C, E, F, G, H, J], circular=True)
    p4 = msu.Path([A, B, ~D, ~C, E, F, G, ~H, J], circular=True)

    paths = {1: p1, 2: p2, 3: p3, 4: p4}
    nodes = {n.id: n for n in [A, B, C, D, E, F, G, H, J]}

    return paths, nodes


def test_find_mergers(generate_core_paths):
    paths, nodes = generate_core_paths

    mg = msu.find_mergers(paths)
    mg_groups = defaultdict(set)
    for source, sink in mg.items():
        mg_groups[sink].add(source)
    sources = list(mg_groups.values())
    assert len(sources) == 4
    assert {"A", "B", "J"} in sources
    assert {"C", "D"} in sources
    assert {"E", "F", "G"} in sources
    assert {"H"} in sources


@pytest.fixture
def load_graph():
    return pp.Pangraph.load_json("tests/data/graph_circ.json")


def test_msu(load_graph):
    pan = load_graph
    MSU_mergers, MSU_paths, MSU_len = msu.minimal_synteny_units(
        pan, L_thr=50, rotate=True
    )
