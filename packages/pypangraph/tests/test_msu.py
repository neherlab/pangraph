import pypangraph as pp
from collections import defaultdict
import pytest

from pypangraph.minimal_synteny_units import (
    minimal_synteny_units,
    flip_msu_to_most_common_orientation,
)
import pypangraph.topology_utils as tu


class TestNode:
    def test_eq(self):
        n1 = tu.OrientedBlock("A", True)
        n2 = tu.OrientedBlock("A", True)
        n3 = tu.OrientedBlock("A", False)
        n4 = tu.OrientedBlock("B", True)
        assert n1 == n2
        assert n1 != n3
        assert n1 != n4

    def test_invert(self):
        n = tu.OrientedBlock("A", True)
        assert n.invert() == tu.OrientedBlock("A", False)


class TestPath:
    def test_constructor_default_nodes_are_not_shared(self):
        p1 = tu.Walk()
        p2 = tu.Walk()
        p1.add_right(tu.OrientedBlock("A", True))
        assert p1.nodes == [tu.OrientedBlock("A", True)]
        assert p2.nodes == []

    def test_eq(self):
        n1 = tu.OrientedBlock("A", True)
        n2 = tu.OrientedBlock("B", True)
        n3 = tu.OrientedBlock("C", True)
        p1 = tu.Walk([n1, n2, n3], circular=True)
        p2 = tu.Walk([n1, n2, n3], circular=True)
        p3 = tu.Walk([n1, ~n2, n3], circular=True)
        p4 = tu.Walk([n1, n2], circular=True)
        assert p1 == p2
        assert p1 != p3
        assert p1 != p4

    def test_invert(self):
        n1 = tu.OrientedBlock("A", True)
        n2 = tu.OrientedBlock("B", True)
        n3 = tu.OrientedBlock("C", True)
        p = tu.Walk([n1, n2, n3], circular=True)
        assert ~p == tu.Walk([~n3, ~n2, ~n1], circular=True)

    def test_rotate_to(self):
        n1 = tu.OrientedBlock("A", True)
        n2 = tu.OrientedBlock("B", True)
        n3 = tu.OrientedBlock("C", True)
        n4 = tu.OrientedBlock("D", True)
        p = tu.Walk([n1, n2, n3, n4], circular=True)
        assert p.rotate_to("B", True) == tu.Walk([n2, n3, n4, n1], circular=True)
        assert p.rotate_to("D", True) == tu.Walk([n4, n1, n2, n3], circular=True)
        assert p.rotate_to("B", False) == tu.Walk([~n2, ~n1, ~n4, ~n3], circular=True)
        assert p.rotate_to("D", False) == tu.Walk([~n4, ~n3, ~n2, ~n1], circular=True)

    def test_rename_bids(self):
        n1 = tu.OrientedBlock("A", True)
        n2 = tu.OrientedBlock("B", False)
        n3 = tu.OrientedBlock("C", True)
        p = tu.Walk([n1, n2, n3], circular=True)
        assert p.rename_bids({"A": "X", "B": "Y", "C": "Z"}) == tu.Walk(
            [tu.OrientedBlock("X", True), tu.OrientedBlock("Y", False), tu.OrientedBlock("Z", True)],
            circular=True,
        )


class TestEdge:
    def test_invert(self):
        n1 = tu.OrientedBlock("A", True)
        n2 = tu.OrientedBlock("B", True)
        e = tu.Edge(n1, n2)
        assert e.invert() == tu.Edge(n2.invert(), n1.invert())

    def test_eq(self):
        n1 = tu.OrientedBlock("A", True)
        n2 = tu.OrientedBlock("B", True)
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
    assert tu.OrientedBlock("A", True).id == "A"
    assert callable(pp.minimal_synteny_units)


@pytest.fixture
def generate_core_paths():
    A = tu.OrientedBlock("A", True)
    B = tu.OrientedBlock("B", True)
    C = tu.OrientedBlock("C", True)  # invert
    D = tu.OrientedBlock("D", True)  # invert
    E = tu.OrientedBlock("E", True)
    F = tu.OrientedBlock("F", True)
    G = tu.OrientedBlock("G", True)
    H = tu.OrientedBlock("H", True)  # invert
    J = tu.OrientedBlock("J", True)

    p1 = tu.Walk([A, B, C, D, E, F, G, H, J], circular=True)
    p2 = tu.Walk([A, B, C, D, E, F, G, H, J], circular=True)
    p3 = tu.Walk([A, B, ~D, ~C, E, F, G, H, J], circular=True)
    p4 = tu.Walk([A, B, ~D, ~C, E, F, G, ~H, J], circular=True)

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
    """Smoke test on the real plasmids graph with invariant checks."""
    pan = load_graph
    MSU_mergers, MSU_paths, MSU_len = minimal_synteny_units(pan, L_thr=50, rotate=True)

    # one MSU walk per strain
    assert set(MSU_paths.keys()) == set(pan.strains())
    # every block maps to a real MSU label
    assert all(label in MSU_len for label in MSU_mergers.values())
    assert all(label.startswith("MSU_") for label in MSU_mergers.values())
    # MSU_0 is the longest unit (labels are assigned by descending length)
    assert max(MSU_len, key=MSU_len.get) == "MSU_0"


def test_minimal_synteny_units_circular(junction_pangraph):
    """End-to-end MSU extraction on the circular junction fixture.

    Core blocks: C1=100(1000bp), C2=200(800bp), C3=300(600bp), C4=400(700bp);
    accessory blocks are dropped by the core filter. Core orders are
    s1/s2 = C1 C2 C3 C4 and s3 = C1 C3 C2 C4. The only edge shared by all three
    strains is the circular wrap C4->C1, so C1 and C4 merge into a single MSU
    (sink=400=C4); C2 and C3 stay singletons. MSUs are labelled by descending
    length: MSU_0={C1,C4} (1700bp), MSU_1=C2 (800bp), MSU_2=C3 (600bp). All paths
    are rotated to MSU_0 (forward); everything is forward so the flip is a no-op,
    and s3 keeps its C2/C3 rearrangement.
    """
    MSU_mergers, MSU_paths, MSU_len = minimal_synteny_units(
        junction_pangraph, L_thr=500, rotate=True
    )

    assert MSU_len == {"MSU_0": 1700, "MSU_1": 800, "MSU_2": 600}
    assert MSU_mergers == {
        "100": "MSU_0",
        "400": "MSU_0",
        "200": "MSU_1",
        "300": "MSU_2",
    }

    expected_s1 = tu.Walk(
        [
            tu.OrientedBlock("MSU_0", True),
            tu.OrientedBlock("MSU_1", True),
            tu.OrientedBlock("MSU_2", True),
        ],
        circular=True,
    )
    expected_s3 = tu.Walk(
        [
            tu.OrientedBlock("MSU_0", True),
            tu.OrientedBlock("MSU_2", True),
            tu.OrientedBlock("MSU_1", True),
        ],
        circular=True,
    )
    assert MSU_paths["s1"] == expected_s1
    assert MSU_paths["s2"] == expected_s1
    assert MSU_paths["s3"] == expected_s3
    assert all(p.circular for p in MSU_paths.values())


def test_minimal_synteny_units_no_rotate(linear_pangraph):
    """MSU extraction with rotate=False on linear paths.

    Core blocks C1=100(1000), C2=200(800), C3=300(600); both strains' core order is
    C1 C2 C3. They form a single transitive chain, so all three merge into one MSU
    (sink=100, 2400bp). Each path reduces to a single MSU node, unrotated.
    """
    MSU_mergers, MSU_paths, MSU_len = minimal_synteny_units(
        linear_pangraph, L_thr=500, rotate=False
    )

    assert MSU_len == {"MSU_0": 2400}
    assert MSU_mergers == {"100": "MSU_0", "200": "MSU_0", "300": "MSU_0"}

    expected = tu.Walk([tu.OrientedBlock("MSU_0", True)], circular=False)
    assert MSU_paths["s1"] == expected
    assert MSU_paths["s2"] == expected
    assert all(p.circular is False for p in MSU_paths.values())


def test_minimal_synteny_units_rotate_requires_circular(linear_pangraph):
    """rotate=True (the default) on linear paths raises a ValueError."""
    with pytest.raises(ValueError, match="Only circular paths"):
        minimal_synteny_units(linear_pangraph, L_thr=500)


def test_flip_msu_to_most_common_orientation():
    """Blocks predominantly on the reverse strand are flipped to forward.

    X is reverse in 2 of 3 walks (net negative) so every X occurrence is flipped;
    Y is forward throughout (net positive) and left untouched. The dict passed in is
    mutated in place and returned.
    """
    paths = {
        "a": tu.Walk(
            [tu.OrientedBlock("X", False), tu.OrientedBlock("Y", True)], circular=True
        ),
        "b": tu.Walk(
            [tu.OrientedBlock("X", False), tu.OrientedBlock("Y", True)], circular=True
        ),
        "c": tu.Walk(
            [tu.OrientedBlock("X", True), tu.OrientedBlock("Y", True)], circular=True
        ),
    }

    result = flip_msu_to_most_common_orientation(paths)

    assert result is paths  # mutates in place and returns the same dict
    # X net = (-1) + (-1) + (+1) = -1 < 0 -> every X occurrence flipped
    assert paths["a"].nodes[0] == tu.OrientedBlock("X", True)
    assert paths["b"].nodes[0] == tu.OrientedBlock("X", True)
    assert paths["c"].nodes[0] == tu.OrientedBlock("X", False)
    # Y net positive -> unchanged
    assert paths["a"].nodes[1] == tu.OrientedBlock("Y", True)
    assert paths["c"].nodes[1] == tu.OrientedBlock("Y", True)
