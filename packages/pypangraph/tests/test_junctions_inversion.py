"""Cross-cutting junction tests on the inversion_pangraph fixture: a circular graph
with a reverse-complemented strain and a single inversion, exercising edges,
positions, stats, and sequences together."""

from Bio.Seq import Seq

import pypangraph.topology_utils as tu
from pypangraph.junctions import BackboneJunctions


def _edge(a, sa, b, sb):
    """Canonical edge string id for OrientedBlock(a, sa) <-> OrientedBlock(b, sb)."""
    return tu.Edge(tu.OrientedBlock(a, sa), tu.OrientedBlock(b, sb)).to_str_id()


def test_inversion_edges_rc_and_private(inversion_pangraph):
    """s1 and its whole-RC s2 share the same 7 backbone edges; s3's inversion adds 2 private ones."""
    bj = BackboneJunctions(inversion_pangraph, L_thr=10)
    freq = bj.stats()["n_isolates"].to_dict()

    shared = {
        _edge("10", True, "50", True),  # C1-C5
        _edge("50", True, "20", True),  # C5-C2
        _edge("20", True, "60", True),  # C2-C6
        _edge("60", True, "30", True),  # C6-C3
        _edge("30", True, "70", True),  # C3-C7
        _edge("70", True, "40", True),  # C7-C4
        _edge("40", True, "10", True),  # C4-C1 (wrap)
    }
    private = {
        _edge("60", True, "70", False),  # C6 -> C7-   (s3 inversion boundary)
        _edge("30", False, "40", True),  # C3- -> C4   (s3 inversion boundary)
    }
    assert set(bj.edges()) == shared | private

    # 5 universal edges (freq 3), 2 broken by the inversion (freq 2), 2 private (freq 1)
    assert sorted(freq.values(), reverse=True) == [3, 3, 3, 3, 3, 2, 2, 1, 1]
    assert all(freq[e] == 1 for e in private)
    assert freq[_edge("60", True, "30", True)] == 2  # C6-C3 broken by inversion
    assert freq[_edge("70", True, "40", True)] == 2  # C7-C4 broken by inversion


def test_inversion_positions_strand_mix(inversion_pangraph):
    """A genome and its reverse complement disagree on canonical strand for every shared edge.

    Unlike junction_pangraph (all core blocks forward, `strand` almost always True), s2 here is
    the whole-genome RC of s1, so for every edge they share the is_canonical flag is opposite.
    """
    bj = BackboneJunctions(inversion_pangraph, L_thr=10)
    pos = bj.positions()

    shared = set(pos.xs("s1", level="iso").index) & set(pos.xs("s2", level="iso").index)
    assert len(shared) == 7
    for edge in shared:
        assert pos.loc[(edge, "s1"), "strand"] != pos.loc[(edge, "s2"), "strand"]

    # genuine mix of both orientations (not the near-constant column of junction_pangraph)
    assert set(pos["strand"]) == {True, False}

    # the inversion flips the canonical sense of the C3-C7 edge in s3 relative to the reference
    c3c7 = _edge("30", True, "70", True)
    assert pos.loc[(c3c7, "s1"), "strand"] != pos.loc[(c3c7, "s3"), "strand"]


def test_inversion_stats(inversion_pangraph):
    """Accessory-bearing edges: a variable center (C5-C2) vs a constant one (C7-C4)."""
    bj = BackboneJunctions(inversion_pangraph, L_thr=10)
    sdf = bj.stats()

    # C5-C2: center [A1] in s1/s2, empty in s3 -> 2 categories, singleton
    row = sdf.loc[_edge("50", True, "20", True)]
    assert row["n_isolates"] == 3
    assert row["n_categories"] == 2
    assert row["n_majority_category"] == 2
    assert not row["is_transitive"]
    assert row["is_singleton"]
    assert row["accessory_length"] == 6  # len(A1)

    # C7-C4: center [A2] in both strains that carry it -> 1 category, transitive
    row = sdf.loc[_edge("70", True, "40", True)]
    assert row["n_isolates"] == 2
    assert row["n_categories"] == 1
    assert row["is_transitive"]
    assert row["accessory_length"] == 8  # len(A2)


def test_inversion_sequences_rc_equivalence(inversion_pangraph):
    """Co-oriented junction sequences of a genome and its RC are identical; reverse-strand center
    blocks are reverse-complemented."""
    bj = BackboneJunctions(inversion_pangraph, L_thr=10)
    # A1 (block 80) consensus, read back from the graph (the block has no per-node edits)
    a1 = next(iter(inversion_pangraph.blocks["80"].to_sequences().values()))

    # C5-C2 carries accessory A1 (reverse in the co-oriented frame) in s1/s2, empty in s3
    seqs = {r.id: str(r.seq) for r in bj.sequences(_edge("50", True, "20", True))}
    assert set(seqs) == {"s1", "s2", "s3"}
    assert seqs["s1"] == seqs["s2"]  # a genome equals its RC once co-oriented
    assert len(seqs["s1"]) == 13 + 6 + 14  # C5 + A1 + C2
    assert len(seqs["s3"]) == 13 + 14  # empty center
    assert (
        str(Seq(a1).reverse_complement()) in seqs["s1"]
    )  # center block reverse-complemented
    assert a1 not in seqs["s1"]

    # C7-C4 (n_isolates=2, s1/s2 only) also co-orients to identical sequences
    seqs2 = {r.id: str(r.seq) for r in bj.sequences(_edge("70", True, "40", True))}
    assert set(seqs2) == {"s1", "s2"}
    assert seqs2["s1"] == seqs2["s2"]
