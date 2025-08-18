import pypangraph as pp
import pypangraph.topology as tp
from pypangraph.msu import minimal_synteny_units
from collections import defaultdict
import pytest


class TestNode:
    def test_eq(self):
        n1 = tp.Node("A", True)
        n2 = tp.Node("A", True)
        n3 = tp.Node("A", False)
        n4 = tp.Node("B", True)
        assert n1 == n2
        assert n1 != n3
        assert n1 != n4

    def test_invert(self):
        n = tp.Node("A", True)
        assert n.invert() == tp.Node("A", False)


class TestPath:
    def test_eq(self):
        n1 = tp.Node("A", True)
        n2 = tp.Node("B", True)
        n3 = tp.Node("C", True)
        p1 = tp.Path([n1, n2, n3], circular=True)
        p2 = tp.Path([n1, n2, n3], circular=True)
        p3 = tp.Path([n1, ~n2, n3], circular=True)
        p4 = tp.Path([n1, n2], circular=True)
        assert p1 == p2
        assert p1 != p3
        assert p1 != p4

    def test_invert(self):
        n1 = tp.Node("A", True)
        n2 = tp.Node("B", True)
        n3 = tp.Node("C", True)
        p = tp.Path([n1, n2, n3], circular=True)
        assert ~p == tp.Path([~n3, ~n2, ~n1], circular=True)

    def test_rotate_to(self):
        n1 = tp.Node("A", True)
        n2 = tp.Node("B", True)
        n3 = tp.Node("C", True)
        n4 = tp.Node("D", True)
        p = tp.Path([n1, n2, n3, n4], circular=True)
        assert p.rotate_to("B", True) == tp.Path([n2, n3, n4, n1], circular=True)
        assert p.rotate_to("D", True) == tp.Path([n4, n1, n2, n3], circular=True)
        assert p.rotate_to("B", False) == tp.Path([~n2, ~n1, ~n4, ~n3], circular=True)
        assert p.rotate_to("D", False) == tp.Path([~n4, ~n3, ~n2, ~n1], circular=True)

    def test_rename_bids(self):
        n1 = tp.Node("A", True)
        n2 = tp.Node("B", False)
        n3 = tp.Node("C", True)
        p = tp.Path([n1, n2, n3], circular=True)
        assert p.rename_bids({"A": "X", "B": "Y", "C": "Z"}) == tp.Path(
            [tp.Node("X", True), tp.Node("Y", False), tp.Node("Z", True)],
            circular=True,
        )


class TestEdge:
    def test_invert(self):
        n1 = tp.Node("A", True)
        n2 = tp.Node("B", True)
        e = tp.Edge(n1, n2)
        assert e.invert() == tp.Edge(n2.invert(), n1.invert())

    def test_eq(self):
        n1 = tp.Node("A", True)
        n2 = tp.Node("B", True)
        e1 = tp.Edge(n1, n2)
        e2 = tp.Edge(n1, n2)
        e3 = tp.Edge(~n1, ~n2)
        e4 = tp.Edge(n1, ~n2)
        e5 = tp.Edge(~n2, ~n1)
        assert e1 == e2
        assert e1 != e3
        assert e1 != e4
        assert e1 == e5


@pytest.fixture
def generate_core_paths():
    A = tp.Node("A", True)
    B = tp.Node("B", True)
    C = tp.Node("C", True)  # invert
    D = tp.Node("D", True)  # invert
    E = tp.Node("E", True)
    F = tp.Node("F", True)
    G = tp.Node("G", True)
    H = tp.Node("H", True)  # invert
    J = tp.Node("J", True)

    p1 = tp.Path([A, B, C, D, E, F, G, H, J], circular=True)
    p2 = tp.Path([A, B, C, D, E, F, G, H, J], circular=True)
    p3 = tp.Path([A, B, ~D, ~C, E, F, G, H, J], circular=True)
    p4 = tp.Path([A, B, ~D, ~C, E, F, G, ~H, J], circular=True)

    paths = {1: p1, 2: p2, 3: p3, 4: p4}
    nodes = {n.id: n for n in [A, B, C, D, E, F, G, H, J]}

    return paths, nodes


def test_find_mergers(generate_core_paths):
    paths, nodes = generate_core_paths

    mg = tp.find_mergers(paths)
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
    return pp.Pangraph.from_json("tests/data/plasmids.json")


def test_msu(load_graph):
    pan = load_graph
    MSU_mergers, MSU_paths, MSU_len = minimal_synteny_units(pan, L_thr=50, rotate=True)
