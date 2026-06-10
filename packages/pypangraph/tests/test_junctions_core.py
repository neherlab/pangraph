"""Tests for the junctions data structures: path splitting, Junction/JunctionNode
primitives, canonical orientation, and the BackboneJunctions interface."""

import pytest

import pypangraph.topology_utils as tu
from pypangraph.junctions import BackboneJunctions
from pypangraph.junctions.junction import Junction, JunctionNode, path_junction_split


# --- path_junction_split tests ---


def test_path_junction_split(junction_pangraph):
    """Splitting s1 at core blocks produces 4 junctions with correct flanking edges."""
    pan = junction_pangraph
    bdf = pan.to_blockstats_df()
    walks = tu.pangraph_to_walks(pan)

    def is_core(bid):
        return (bdf.loc[bid, "len"] >= 500) and bdf.loc[bid, "core"]

    junctions = path_junction_split(walks["s1"], is_core)
    assert len(junctions) == 4

    edges = [j.flanking_edge().to_str_id() for j in junctions]
    assert set(edges) == {
        "100_r__400_r",
        "100_f__200_f",
        "200_f__300_f",
        "300_f__400_f",
    }

    # The junction between C1 and C2 should contain 2 accessory blocks (A1, A2)
    j_c1_c2 = [j for j in junctions if j.flanking_edge().to_str_id() == "100_f__200_f"][
        0
    ]
    assert len(j_c1_c2.center) == 2


def test_path_junction_split_rearranged(junction_pangraph):
    """Splitting s3 (rearranged) produces 4 junctions with different edges than s1/s2."""
    pan = junction_pangraph
    bdf = pan.to_blockstats_df()
    walks = tu.pangraph_to_walks(pan)

    def is_core(bid):
        return (bdf.loc[bid, "len"] >= 500) and bdf.loc[bid, "core"]

    junctions = path_junction_split(walks["s3"], is_core)
    assert len(junctions) == 4

    edges = {j.flanking_edge().to_str_id() for j in junctions}
    assert edges == {"100_r__400_r", "100_f__300_f", "200_r__300_r", "200_f__400_f"}


def test_path_junction_split_requires_two_core_blocks():
    """path_junction_split errors when a path has fewer than 2 core blocks.

    A circular path with no core blocks previously crashed with IndexError; a path
    with a single core block cannot form a junction either.
    """
    # no core blocks (circular) — would otherwise hit IndexError on junctions[0]
    p0 = tu.Walk([tu.OrientedBlock(1, True), tu.OrientedBlock(2, True)], circular=True)
    with pytest.raises(ValueError, match="at least 2"):
        path_junction_split(p0, lambda bid: False)

    # exactly one core block is still insufficient
    p1 = tu.Walk([tu.OrientedBlock(1, True), tu.OrientedBlock(2, True)], circular=True)
    with pytest.raises(ValueError, match="at least 2"):
        path_junction_split(p1, lambda bid: bid == 1)


def test_path_junction_split_linear(linear_pangraph):
    """Splitting linear paths produces terminal junctions with None flanks.

    s1: c5+ C1+ A2+ C2+ C3+
    s2: C1+ A3+ C2+ C3+ c5+
    """
    pan = linear_pangraph
    bdf = pan.to_blockstats_df()
    walks = tu.pangraph_to_walks(pan)

    def is_core(bid):
        return (bdf.loc[bid, "len"] >= 500) and bdf.loc[bid, "core"]

    # s1: c5+ C1+ A2+ C2+ C3+
    # Expected: [None|c5+|C1+], [C1+||C2+], [C2+||C3+], no trailing (C3 is last)
    # Wait - trailing after C3 is empty, but left_node=C3 so we get [C3+||None]
    junctions_s1 = path_junction_split(walks["s1"], is_core)

    # Leading terminal: left=None, center=[c5+], right=C1+
    assert junctions_s1[0].left is None
    assert len(junctions_s1[0].center) == 1
    assert junctions_s1[0].right is not None

    # Internal junctions have both flanks
    assert junctions_s1[1].left is not None
    assert junctions_s1[1].right is not None

    # Trailing terminal: left=C3+, center=[], right=None
    assert junctions_s1[-1].left is not None
    assert junctions_s1[-1].right is None

    # flanking_edge is None for terminal junctions
    assert junctions_s1[0].flanking_edge() is None
    assert junctions_s1[-1].flanking_edge() is None
    assert junctions_s1[1].flanking_edge() is not None

    # s2: C1+ A3+ C2+ C3+ c5+
    junctions_s2 = path_junction_split(walks["s2"], is_core)

    # No leading terminal (C1 is first)
    assert junctions_s2[0].left is None
    assert len(junctions_s2[0].center) == 0
    assert junctions_s2[0].right is not None

    # Trailing terminal: left=C3+, center=[c5+], right=None
    assert junctions_s2[-1].left is not None
    assert len(junctions_s2[-1].center) == 1
    assert junctions_s2[-1].right is None


# --- JunctionNode / Junction primitives ---


def test_junction_node_inherits_equality():
    """JunctionNode equality uses block_id + strand only (ignores node_id)."""
    a = JunctionNode(100, True, 1)
    b = JunctionNode(100, True, 2)
    c = JunctionNode(100, False, 1)
    assert a == b
    assert a != c
    assert hash(a) == hash(b)


def test_junction_node_invert():
    """JunctionNode.invert() flips strand but preserves node_id."""
    n = JunctionNode(100, True, 42)
    inv = n.invert()
    assert inv.id == 100
    assert inv.strand is False
    assert inv.node_id == 42
    assert isinstance(inv, JunctionNode)


def test_junction_oriented_blocks():
    """Junction.oriented_blocks() flattens left + center + right, skipping None flanks."""
    left = tu.OrientedBlock(100, True)
    a1 = tu.OrientedBlock(200, True)
    a2 = tu.OrientedBlock(300, False)
    right = tu.OrientedBlock(400, True)

    # full junction: left + non-empty center + right
    full = Junction(left, tu.Walk([a1, a2]), right)
    assert full.oriented_blocks() == [left, a1, a2, right]

    # empty center: just the two flanks
    empty_center = Junction(left, tu.Walk([]), right)
    assert empty_center.oriented_blocks() == [left, right]

    # terminal junction with no left flank
    no_left = Junction(None, tu.Walk([a1, a2]), right)
    assert no_left.oriented_blocks() == [a1, a2, right]

    # terminal junction with no right flank
    no_right = Junction(left, tu.Walk([a1, a2]), None)
    assert no_right.oriented_blocks() == [left, a1, a2]


def test_junction_invert_terminal():
    """invert() is robust to terminal junctions: a None flank stays None and the
    present flank moves to the opposite side with its node_id preserved."""
    a1 = tu.OrientedBlock(200, True)
    a2 = tu.OrientedBlock(300, False)

    # no left flank -> after inversion, no right flank
    right = JunctionNode(400, True, 42)
    no_left = Junction(None, tu.Walk([a1, a2]), right)
    inv = no_left.invert()
    assert inv.right is None
    assert inv.left == right.invert()
    assert inv.left.node_id == 42  # node_id preserved across inversion
    assert inv.center == tu.Walk([a2.invert(), a1.invert()])
    # double inversion round-trips
    assert no_left.invert().invert() == no_left

    # symmetric: no right flank -> after inversion, no left flank
    left = JunctionNode(100, True, 7)
    no_right = Junction(left, tu.Walk([a1, a2]), None)
    rinv = no_right.invert()
    assert rinv.left is None
    assert rinv.right == left.invert()
    assert rinv.right.node_id == 7


# --- is_canonical / to_canonical ---


def test_edge_is_canonical():
    """Edge.is_canonical() is True iff (left, right) is the lex-min orientation."""
    # natural form "100_f__200_f" is already lex-min
    e = tu.Edge.from_str_id("100_f__200_f")
    assert e.is_canonical()
    assert e.to_str_id() == "100_f__200_f"

    # inverting flips the orientation: now non-canonical
    inv = e.invert()
    assert not inv.is_canonical()
    # but to_str_id() always returns the canonical form
    assert inv.to_str_id() == "100_f__200_f"

    # RC-palindromic edge: invert() yields the same string; tie resolves canonical
    palindrome = tu.Edge(tu.OrientedBlock(100, True), tu.OrientedBlock(100, False))
    assert palindrome.is_canonical()


def test_str_id_roundtrip_underscore_ids():
    """to_str_id / from_str_id round-trip block ids that themselves contain '_'
    (e.g. MSU-renamed ids like 'MSU_0'); the strand suffix is peeled from the right."""
    ob = tu.OrientedBlock("MSU_0", True)
    assert tu.OrientedBlock.from_str_id(ob.to_str_id()) == ob

    e = tu.Edge(tu.OrientedBlock("MSU_0", True), tu.OrientedBlock("MSU_12", False))
    assert tu.Edge.from_str_id(e.to_str_id()) == e


def test_junction_is_canonical(junction_pangraph):
    """Junction.is_canonical() reports whether the junction's natural orientation
    matches the lex-min canonical form of its flanking edge."""
    bj = BackboneJunctions(junction_pangraph, L_thr=500)

    # s1 walks C1+ ... C2+: natural form "100_f__200_f" == canonical
    assert bj["100_f__200_f"]["s1"].is_canonical()

    # s3 walks C3+ then C2+ for the C2-C3 adjacency: natural form
    # "300_f__200_f" loses lex-min to its invert "200_r__300_r"
    assert not bj["200_r__300_r"]["s3"].is_canonical()


def test_junction_to_canonical(junction_pangraph):
    """to_canonical() returns self when already canonical, inverted otherwise;
    the result always satisfies is_canonical()."""
    bj = BackboneJunctions(junction_pangraph, L_thr=500)

    # canonical case: to_canonical() is a no-op
    j_canon = bj["100_f__200_f"]["s1"]
    assert j_canon.to_canonical() is j_canon

    # non-canonical case: to_canonical() inverts and the result is canonical
    j_inv = bj["200_r__300_r"]["s3"]
    j_out = j_inv.to_canonical()
    assert j_out is not j_inv
    assert j_out.is_canonical()
    assert j_out == j_inv.invert()


def test_junction_canonical_terminal_raises():
    """is_canonical()/to_canonical() are undefined for terminal junctions and raise,
    so they cannot be applied accidentally to a junction with a missing flank."""
    center = tu.Walk([tu.OrientedBlock(200, True)])

    no_left = Junction(None, center, tu.OrientedBlock(400, True))
    with pytest.raises(ValueError, match="terminal"):
        no_left.is_canonical()
    with pytest.raises(ValueError, match="terminal"):
        no_left.to_canonical()

    no_right = Junction(tu.OrientedBlock(100, True), center, None)
    with pytest.raises(ValueError, match="terminal"):
        no_right.is_canonical()
    with pytest.raises(ValueError, match="terminal"):
        no_right.to_canonical()


# --- BackboneJunctions interface ---


def test_junctions_edge_freq(junction_pangraph):
    """Edge isolate count is correct and sorted descending."""
    bj = BackboneJunctions(junction_pangraph, L_thr=500)
    edge_freq = bj.stats()["n_isolates"]

    # C4→C1 edge is universal (all 3 strains)
    assert edge_freq["100_r__400_r"] == 3
    # C1→C2, C2→C3, C3→C4 edges shared by s1 and s2
    assert edge_freq["100_f__200_f"] == 2
    assert edge_freq["200_f__300_f"] == 2
    assert edge_freq["300_f__400_f"] == 2
    # s3-private edges
    assert edge_freq["100_f__300_f"] == 1
    assert edge_freq["200_f__400_f"] == 1
    assert edge_freq["200_r__300_r"] == 1

    # sorted descending
    counts = edge_freq.values
    assert all(counts[i] >= counts[i + 1] for i in range(len(counts) - 1))


def test_backbone_getitem(junction_pangraph):
    """bj[edge_str] returns the {isolate -> Junction} mapping for an edge,
    and `edge in bj` reports presence."""
    bj = BackboneJunctions(junction_pangraph, L_thr=500)

    per_iso = bj["100_f__200_f"]
    assert set(per_iso.keys()) == {"s1", "s2"}
    j = per_iso["s1"]
    assert j.flanking_edge().to_str_id() == "100_f__200_f"
    assert len(j.center) == 2  # A1 + A2

    assert "100_f__200_f" in bj
    assert "nonexistent_edge" not in bj
    with pytest.raises(KeyError):
        _ = bj["nonexistent_edge"]


def test_backbone_edges(junction_pangraph):
    """edges() returns all distinct edge string IDs."""
    bj = BackboneJunctions(junction_pangraph, L_thr=500)
    assert set(bj.edges()) == {
        "100_r__400_r",
        "100_f__200_f",
        "200_f__300_f",
        "300_f__400_f",
        "100_f__300_f",
        "200_f__400_f",
        "200_r__300_r",
    }
