import pypangraph as pp
from collections import defaultdict
import pytest

from pypangraph.minimal_synteny_units import minimal_synteny_units
import pypangraph.topology_utils as tu


class TestNode:
    def test_eq(self):
        n1 = tu.Node("A", True)
        n2 = tu.Node("A", True)
        n3 = tu.Node("A", False)
        n4 = tu.Node("B", True)
        assert n1 == n2
        assert n1 != n3
        assert n1 != n4

    def test_invert(self):
        n = tu.Node("A", True)
        assert n.invert() == tu.Node("A", False)


class TestPath:
    def test_constructor_default_nodes_are_not_shared(self):
        p1 = tu.Path()
        p2 = tu.Path()
        p1.add_right(tu.Node("A", True))
        assert p1.nodes == [tu.Node("A", True)]
        assert p2.nodes == []

    def test_eq(self):
        n1 = tu.Node("A", True)
        n2 = tu.Node("B", True)
        n3 = tu.Node("C", True)
        p1 = tu.Path([n1, n2, n3], circular=True)
        p2 = tu.Path([n1, n2, n3], circular=True)
        p3 = tu.Path([n1, ~n2, n3], circular=True)
        p4 = tu.Path([n1, n2], circular=True)
        assert p1 == p2
        assert p1 != p3
        assert p1 != p4

    def test_invert(self):
        n1 = tu.Node("A", True)
        n2 = tu.Node("B", True)
        n3 = tu.Node("C", True)
        p = tu.Path([n1, n2, n3], circular=True)
        assert ~p == tu.Path([~n3, ~n2, ~n1], circular=True)

    def test_rotate_to(self):
        n1 = tu.Node("A", True)
        n2 = tu.Node("B", True)
        n3 = tu.Node("C", True)
        n4 = tu.Node("D", True)
        p = tu.Path([n1, n2, n3, n4], circular=True)
        assert p.rotate_to("B", True) == tu.Path([n2, n3, n4, n1], circular=True)
        assert p.rotate_to("D", True) == tu.Path([n4, n1, n2, n3], circular=True)
        assert p.rotate_to("B", False) == tu.Path([~n2, ~n1, ~n4, ~n3], circular=True)
        assert p.rotate_to("D", False) == tu.Path([~n4, ~n3, ~n2, ~n1], circular=True)

    def test_rename_bids(self):
        n1 = tu.Node("A", True)
        n2 = tu.Node("B", False)
        n3 = tu.Node("C", True)
        p = tu.Path([n1, n2, n3], circular=True)
        assert p.rename_bids({"A": "X", "B": "Y", "C": "Z"}) == tu.Path(
            [tu.Node("X", True), tu.Node("Y", False), tu.Node("Z", True)],
            circular=True,
        )


class TestEdge:
    def test_invert(self):
        n1 = tu.Node("A", True)
        n2 = tu.Node("B", True)
        e = tu.Edge(n1, n2)
        assert e.invert() == tu.Edge(n2.invert(), n1.invert())

    def test_eq(self):
        n1 = tu.Node("A", True)
        n2 = tu.Node("B", True)
        e1 = tu.Edge(n1, n2)
        e2 = tu.Edge(n1, n2)
        e3 = tu.Edge(~n1, ~n2)
        e4 = tu.Edge(n1, ~n2)
        e5 = tu.Edge(~n2, ~n1)
        assert e1 == e2
        assert e1 != e3
        assert e1 != e4
        assert e1 == e5


def test_new_import_paths_smoke():
    assert callable(minimal_synteny_units)
    assert tu.Node("A", True).id == "A"
    assert callable(pp.minimal_synteny_units)


@pytest.fixture
def generate_core_paths():
    A = tu.Node("A", True)
    B = tu.Node("B", True)
    C = tu.Node("C", True)  # invert
    D = tu.Node("D", True)  # invert
    E = tu.Node("E", True)
    F = tu.Node("F", True)
    G = tu.Node("G", True)
    H = tu.Node("H", True)  # invert
    J = tu.Node("J", True)

    p1 = tu.Path([A, B, C, D, E, F, G, H, J], circular=True)
    p2 = tu.Path([A, B, C, D, E, F, G, H, J], circular=True)
    p3 = tu.Path([A, B, ~D, ~C, E, F, G, H, J], circular=True)
    p4 = tu.Path([A, B, ~D, ~C, E, F, G, ~H, J], circular=True)

    paths = {1: p1, 2: p2, 3: p3, 4: p4}
    nodes = {n.id: n for n in [A, B, C, D, E, F, G, H, J]}

    return paths, nodes


def test_find_mergers(generate_core_paths):
    paths, nodes = generate_core_paths

    mg = tu.find_mergers(paths)
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
