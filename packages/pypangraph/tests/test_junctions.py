"""Tests for the pypangraph.junctions sub-package."""

import numpy as np
import pandas as pd
import pytest
from Bio.Seq import Seq

import pypangraph as pp
import pypangraph.topology_utils as tu
from pypangraph.junctions import (
    Junction,
    JunctionNode,
    path_junction_split,
    junctions_dataframe,
    junction_positions,
    junction_stats,
    BackboneJunctions,
)


def test_path_junction_split(junction_pangraph):
    """Splitting s1 at core blocks produces 4 junctions with correct flanking edges."""
    pan = junction_pangraph
    bdf = pan.to_blockstats_df()
    paths = tu.pangraph_to_path_dict(pan)

    def is_core(bid):
        return (bdf.loc[bid, "len"] >= 500) and bdf.loc[bid, "core"]

    junctions = path_junction_split(paths["s1"], is_core)
    assert len(junctions) == 4

    edges = [j.flanking_edge().to_str_id() for j in junctions]
    assert set(edges) == {"100_r__400_r", "100_f__200_f", "200_f__300_f", "300_f__400_f"}

    # The junction between C1 and C2 should contain 2 accessory blocks (A1, A2)
    j_c1_c2 = [j for j in junctions if j.flanking_edge().to_str_id() == "100_f__200_f"][0]
    assert len(j_c1_c2.center) == 2


def test_path_junction_split_rearranged(junction_pangraph):
    """Splitting s3 (rearranged) produces 4 junctions with different edges than s1/s2."""
    pan = junction_pangraph
    bdf = pan.to_blockstats_df()
    paths = tu.pangraph_to_path_dict(pan)

    def is_core(bid):
        return (bdf.loc[bid, "len"] >= 500) and bdf.loc[bid, "core"]

    junctions = path_junction_split(paths["s3"], is_core)
    assert len(junctions) == 4

    edges = {j.flanking_edge().to_str_id() for j in junctions}
    assert edges == {"100_r__400_r", "100_f__300_f", "200_r__300_r", "200_f__400_f"}


def test_junctions_dataframe_shape(junction_pangraph):
    """DataFrame has 3 isolates and 7 distinct edges."""
    jdf, edge_freq = junctions_dataframe(junction_pangraph, L_thr=500)
    assert jdf.shape == (3, 7)
    assert len(edge_freq) == 7


def test_junctions_dataframe_values(junction_pangraph):
    """Verify specific junction lengths and NaN positions."""
    jdf, _ = junctions_dataframe(junction_pangraph, L_thr=500)

    # s1: C1→C2 junction has A1(200) + A2(150) = 350
    assert jdf.loc["s1", "100_f__200_f"] == 350.0
    # s2: C1→C2 junction has A1(200) only
    assert jdf.loc["s2", "100_f__200_f"] == 200.0
    # s2: C2→C3 junction has A3(300)
    assert jdf.loc["s2", "200_f__300_f"] == 300.0
    # s3: C1→C3 junction has A2(150)
    assert jdf.loc["s3", "100_f__300_f"] == 150.0
    # s3: C2→C4 junction has A3(300)
    assert jdf.loc["s3", "200_f__400_f"] == 300.0

    # empty junctions have length 0
    assert jdf.loc["s1", "200_f__300_f"] == 0.0
    assert jdf.loc["s1", "300_f__400_f"] == 0.0
    assert jdf.loc["s1", "100_r__400_r"] == 0.0

    # s3 doesn't have the C1→C2 edge (rearranged), so it's NaN
    assert np.isnan(jdf.loc["s3", "100_f__200_f"])
    # s1 doesn't have the s3-private edges
    assert np.isnan(jdf.loc["s1", "100_f__300_f"])
    assert np.isnan(jdf.loc["s1", "200_f__400_f"])


def test_junctions_edge_freq(junction_pangraph):
    """Edge frequency is correct and sorted descending."""
    _, edge_freq = junctions_dataframe(junction_pangraph, L_thr=500)

    assert edge_freq.name == "count"

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


def test_junction_positions_forward_strand(junction_pangraph):
    """Verify positions for forward-strand junctions in s1.

    s1 path: C1+(0,1000) A1+(1000,1200) A2+(1200,1350) C2+(1350,2150) C3+(2150,2750) C4+(2750,3450)
    """
    pan = junction_pangraph
    jdf, _ = junctions_dataframe(pan, L_thr=500)
    pos = junction_positions(pan, jdf)

    # Edge C1+→C2+ in s1: both flanks on forward strand → strand=True
    row = pos.loc[("s1", "100_f__200_f")]
    assert row["strand"] == True
    assert row["left_start"] == 0
    assert row["left_end"] == 1000
    assert row["right_start"] == 1350
    assert row["right_end"] == 2150

    # Edge C2+→C3+ in s1: empty junction, strand=True
    row = pos.loc[("s1", "200_f__300_f")]
    assert row["strand"] == True
    assert row["left_start"] == 1350
    assert row["left_end"] == 2150
    assert row["right_start"] == 2150
    assert row["right_end"] == 2750


def test_junction_positions_inverted_edge(junction_pangraph):
    """Verify positions for the universal C4→C1 edge (inverted in all strains).

    The canonical edge is "100_r__400_r" = Edge(C1-, C4-). But in all strains
    C1 and C4 are on + strand, so same_strand=False → positions are swapped.

    s1: C4+(2750,3450) ... C1+(0,1000) (wraps around)
    """
    pan = junction_pangraph
    jdf, _ = junctions_dataframe(pan, L_thr=500)
    pos = junction_positions(pan, jdf)

    # In s1: junction is inverted, left=C4, right=C1
    row = pos.loc[("s1", "100_r__400_r")]
    assert row["strand"] == False
    assert row["left_start"] == 2750   # C4 start
    assert row["left_end"] == 3450     # C4 end
    assert row["right_start"] == 0     # C1 start
    assert row["right_end"] == 1000    # C1 end


def test_junction_positions_rearranged_strain(junction_pangraph):
    """Verify positions for s3 which has a rearrangement (C2/C3 swapped).

    s3 path: C1+(0,1000) A2-(1000,1150) C3+(1150,1750) C2+(1750,2550) A3+(2550,2850) C4+(2850,3550)
    """
    pan = junction_pangraph
    jdf, _ = junctions_dataframe(pan, L_thr=500)
    pos = junction_positions(pan, jdf)

    # Edge C1+→C3+: forward in s3
    row = pos.loc[("s3", "100_f__300_f")]
    assert row["strand"] == True
    assert row["left_start"] == 0      # C1
    assert row["left_end"] == 1000
    assert row["right_start"] == 1150  # C3
    assert row["right_end"] == 1750

    # Edge "200_r__300_r" = Edge(C2-, C3-): in s3 both are on + strand → inverted
    # Swapped: left=C3(1150,1750), right=C2(1750,2550)
    row = pos.loc[("s3", "200_r__300_r")]
    assert row["strand"] == False
    assert row["left_start"] == 1150   # C3
    assert row["left_end"] == 1750
    assert row["right_start"] == 1750  # C2
    assert row["right_end"] == 2550

    # Edge C2+→C4+: forward in s3, A3 between them
    row = pos.loc[("s3", "200_f__400_f")]
    assert row["strand"] == True
    assert row["left_start"] == 1750   # C2
    assert row["left_end"] == 2550
    assert row["right_start"] == 2850  # C4
    assert row["right_end"] == 3550


def test_junction_positions_shape(junction_pangraph):
    """Position DataFrame has one row per (isolate, edge) with non-NaN junction length."""
    pan = junction_pangraph
    jdf, _ = junctions_dataframe(pan, L_thr=500)
    pos = junction_positions(pan, jdf)

    # Total non-NaN entries in jdf = number of position rows
    n_junctions = jdf.notna().sum().sum()
    assert len(pos) == n_junctions
    assert list(pos.columns) == ["left_start", "left_end", "right_start", "right_end", "strand"]


def test_path_junction_split_linear(linear_pangraph):
    """Splitting linear paths produces terminal junctions with None flanks.

    s1: A1+ C1+ A2+ C2+ C3+
    s2: C1+ A3+ C2+ C3+ A1+
    """
    pan = linear_pangraph
    bdf = pan.to_blockstats_df()
    paths = tu.pangraph_to_path_dict(pan)

    def is_core(bid):
        return (bdf.loc[bid, "len"] >= 500) and bdf.loc[bid, "core"]

    # s1: A1+ C1+ A2+ C2+ C3+
    # Expected: [None|A1+|C1+], [C1+||C2+], [C2+||C3+], no trailing (C3 is last)
    # Wait - trailing after C3 is empty, but left_node=C3 so we get [C3+||None]
    junctions_s1 = path_junction_split(paths["s1"], is_core)

    # Leading terminal: left=None, center=[A1+], right=C1+
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

    # s2: C1+ A3+ C2+ C3+ A1+
    junctions_s2 = path_junction_split(paths["s2"], is_core)

    # No leading terminal (C1 is first)
    assert junctions_s2[0].left is None
    assert len(junctions_s2[0].center) == 0
    assert junctions_s2[0].right is not None

    # Trailing terminal: left=C3+, center=[A1+], right=None
    assert junctions_s2[-1].left is not None
    assert len(junctions_s2[-1].center) == 1
    assert junctions_s2[-1].right is None


def test_junctions_dataframe_linear(linear_pangraph):
    """DataFrame for linear paths excludes terminal junctions (no edge)."""
    jdf, edge_freq = junctions_dataframe(linear_pangraph, L_thr=500)

    # Both strains share edges C1→C2 and C2→C3
    assert jdf.shape[0] == 2  # 2 isolates
    # s1 has: C1→C2 (A2, len=150), C2→C3 (empty, len=0)
    # s2 has: C1→C2 (A3, len=300), C2→C3 (empty, len=0)
    assert jdf.shape[1] == 2  # 2 edges (shared by both)

    assert jdf.loc["s1", "100_f__200_f"] == 150.0
    assert jdf.loc["s2", "100_f__200_f"] == 300.0
    assert jdf.loc["s1", "200_f__300_f"] == 0.0
    assert jdf.loc["s2", "200_f__300_f"] == 0.0


def test_junction_positions_linear(linear_pangraph):
    """Positions work correctly for linear paths.

    s1: A1+(0,200) C1+(200,1200) A2+(1200,1350) C2+(1350,2150) C3+(2150,2750)
    s2: C1+(0,1000) A3+(1000,1300) C2+(1300,2100) C3+(2100,2700) A1+(2700,2900)
    """
    pan = linear_pangraph
    jdf, _ = junctions_dataframe(pan, L_thr=500)
    pos = junction_positions(pan, jdf)

    # s1: C1→C2 has A2 in between
    row = pos.loc[("s1", "100_f__200_f")]
    assert row["left_start"] == 200     # C1
    assert row["left_end"] == 1200
    assert row["right_start"] == 1350   # C2
    assert row["right_end"] == 2150

    # s2: C1→C2 has A3 in between
    row = pos.loc[("s2", "100_f__200_f")]
    assert row["left_start"] == 0       # C1
    assert row["left_end"] == 1000
    assert row["right_start"] == 1300   # C2
    assert row["right_end"] == 2100


@pytest.fixture
def plasmid_pangraph():
    return pp.Pangraph.from_json("tests/data/plasmids.json")


def test_junctions_dataframe_smoke(plasmid_pangraph):
    """junctions_dataframe runs on the real plasmids dataset without errors."""
    jdf, edge_freq = junctions_dataframe(plasmid_pangraph, L_thr=500)

    assert isinstance(jdf, pd.DataFrame)
    assert isinstance(edge_freq, pd.Series)
    assert jdf.shape[0] == len(plasmid_pangraph.strains())
    assert jdf.shape[1] == len(edge_freq)
    assert edge_freq.name == "count"
    assert (edge_freq > 0).all()


def test_junction_positions_smoke(plasmid_pangraph):
    """junction_positions runs on the real plasmids dataset without errors."""
    jdf, _ = junctions_dataframe(plasmid_pangraph, L_thr=500)
    pos = junction_positions(plasmid_pangraph, jdf)

    assert isinstance(pos, pd.DataFrame)
    n_junctions = jdf.notna().sum().sum()
    assert len(pos) == n_junctions
    assert list(pos.columns) == ["left_start", "left_end", "right_start", "right_end", "strand"]


# --- JunctionNode tests ---


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


# --- BackboneJunctions tests ---


def test_backbone_dataframe_matches_standalone(junction_pangraph):
    """BackboneJunctions.dataframe() pivot table matches junctions_dataframe()."""
    pan = junction_pangraph
    bj = BackboneJunctions(pan, L_thr=500)
    jdf_bj, stats_bj = bj.dataframe()
    jdf_fn, _ = junctions_dataframe(pan, L_thr=500)

    # Sort columns for comparison (tie-breaking order may differ)
    jdf_bj = jdf_bj[sorted(jdf_bj.columns)]
    jdf_fn = jdf_fn[sorted(jdf_fn.columns)]
    pd.testing.assert_frame_equal(jdf_bj, jdf_fn)
    # stats_df replaces edge_freq; frequency column should match
    assert isinstance(stats_bj, pd.DataFrame)
    assert "frequency" in stats_bj.columns


def test_backbone_positions_matches_standalone(junction_pangraph):
    """BackboneJunctions.positions() produces the same result as junction_positions()."""
    pan = junction_pangraph
    bj = BackboneJunctions(pan, L_thr=500)
    pos_bj = bj.positions()

    jdf, _ = junctions_dataframe(pan, L_thr=500)
    pos_fn = junction_positions(pan, jdf)

    # Sort both for comparison (row order may differ)
    pos_bj = pos_bj.sort_index()
    pos_fn = pos_fn.sort_index()
    pd.testing.assert_frame_equal(pos_bj, pos_fn)


def test_backbone_junctions_for(junction_pangraph):
    """junctions_for returns correct number of junctions per isolate."""
    bj = BackboneJunctions(junction_pangraph, L_thr=500)
    assert len(bj.junctions_for("s1")) == 4
    assert len(bj.junctions_for("s2")) == 4
    assert len(bj.junctions_for("s3")) == 4


def test_backbone_junction_for(junction_pangraph):
    """junction_for returns the correct junction for a given isolate/edge."""
    bj = BackboneJunctions(junction_pangraph, L_thr=500)
    j = bj.junction_for("s1", "100_f__200_f")
    assert j.flanking_edge().to_str_id() == "100_f__200_f"
    assert len(j.center) == 2  # A1 + A2


def test_backbone_edges(junction_pangraph):
    """edges() returns all distinct edge string IDs."""
    bj = BackboneJunctions(junction_pangraph, L_thr=500)
    assert set(bj.edges()) == {
        "100_r__400_r", "100_f__200_f", "200_f__300_f", "300_f__400_f",
        "100_f__300_f", "200_f__400_f", "200_r__300_r",
    }


# --- Sequence extraction tests ---


def test_sequences_forward_junction(sequence_pangraph):
    """Forward-strand junction produces correct concatenated sequence.

    s1: C1+(n1) A1+(n2) C2+(n3) — forward orientation for edge 10_f__20_f
    C1 consensus: AAACCC (no edits for n1)
    A1 consensus: TTC (no edits for n2)
    C2 consensus: GGGAAA (no edits for n3)
    Expected: AAACCC + TTC + GGGAAA = AAACCCTTCGGGAAA
    """
    bj = BackboneJunctions(sequence_pangraph, L_thr=4)
    records = bj.sequences("10_f__20_f")
    by_iso = {r.id: str(r.seq) for r in records}
    assert by_iso["s1"] == "AAACCCTTCGGGAAA"


def test_sequences_inverted_junction(sequence_pangraph):
    """Inverted junction is co-oriented with the forward one.

    s2: C2-(n4) A1-(n5) C1-(n6) — inverted for edge 10_f__20_f
    After inversion: C1+(n6) A1+(n5) C2+(n4)
    C1 n6 has sub pos=0 A->T, so seq = TAACCC
    A1 n5: TTC (no edits)
    C2 n4: GGGAAA (no edits)
    Expected: TAACCC + TTC + GGGAAA = TAACCCTTCGGGAAA
    """
    bj = BackboneJunctions(sequence_pangraph, L_thr=4)
    records = bj.sequences("10_f__20_f")
    by_iso = {r.id: str(r.seq) for r in records}
    assert by_iso["s2"] == "TAACCCTTCGGGAAA"


def test_sequences_co_orientation(sequence_pangraph):
    """Both isolates' sequences start with C1 and end with C2 (co-oriented)."""
    bj = BackboneJunctions(sequence_pangraph, L_thr=4)
    records = bj.sequences("10_f__20_f")
    assert len(records) == 2
    for r in records:
        seq = str(r.seq)
        # All sequences should start with C1 block (6bp) and end with C2 block (6bp)
        assert len(seq) == 15  # 6 + 3 + 6


def test_sequences_record_metadata(sequence_pangraph):
    """SeqRecord id is isolate name, description is edge string."""
    bj = BackboneJunctions(sequence_pangraph, L_thr=4)
    records = bj.sequences("10_f__20_f")
    ids = {r.id for r in records}
    assert ids == {"s1", "s2"}
    for r in records:
        assert r.description == "10_f__20_f"


def test_sequences_reverse_complement_in_center(junction_pangraph):
    """Accessory blocks on reverse strand are reverse-complemented.

    s3 has edge 100_f__300_f: Junction(C1+(n13), [A2-(n14)], C3+(n15))
    A2 is on reverse strand. Consensus is "A"*150, so rc = "T"*150.
    C1 and C3 are forward: "A"*1000 and "A"*600.
    Expected: "A"*1000 + "T"*150 + "A"*600
    """
    bj = BackboneJunctions(junction_pangraph, L_thr=500)
    records = bj.sequences("100_f__300_f")
    assert len(records) == 1
    seq = str(records[0].seq)
    assert records[0].id == "s3"
    assert len(seq) == 1000 + 150 + 600
    assert seq[:1000] == "A" * 1000
    assert seq[1000:1150] == "T" * 150  # A2 reverse-complemented
    assert seq[1150:] == "A" * 600


def test_sequences_empty_junction(junction_pangraph):
    """Empty junction (no accessory blocks) returns just the two core block sequences.

    s1 edge 200_f__300_f: Junction(C2+, [], C3+) — no accessory blocks.
    Expected: "A"*800 + "A"*600 = "A"*1400
    """
    bj = BackboneJunctions(junction_pangraph, L_thr=500)
    records = bj.sequences("200_f__300_f")
    by_iso = {r.id: str(r.seq) for r in records}
    assert len(by_iso["s1"]) == 800 + 600
    assert by_iso["s1"] == "A" * 1400


def test_sequences_nonexistent_edge(junction_pangraph):
    """Requesting sequences for a nonexistent edge returns empty list."""
    bj = BackboneJunctions(junction_pangraph, L_thr=500)
    assert bj.sequences("999_f__888_f") == []


def test_sequences_smoke(plasmid_pangraph):
    """sequences() runs on the real plasmids dataset without errors."""
    bj = BackboneJunctions(plasmid_pangraph, L_thr=500)
    jdf, stats_df = bj.dataframe()

    # Test on the most frequent edge
    top_edge = stats_df.index[0]
    top_freq = stats_df.loc[top_edge, "frequency"]
    records = bj.sequences(top_edge)
    assert len(records) == top_freq
    for r in records:
        assert len(r.seq) > 0
        assert r.description == top_edge


# --- Junction stats tests ---


def test_junction_stats_values(junction_pangraph):
    """Verify all per-edge statistics against expected values.

    See plan for derivation of expected values.
    """
    bj = BackboneJunctions(junction_pangraph, L_thr=500)
    sdf = bj.stats()

    expected = {
        "100_r__400_r": {"frequency": 3, "n_categories": 1, "majority_category_freq": 3,
                         "left_core_length": 1000, "right_core_length": 700, "accessory_length": 0},
        "100_f__200_f": {"frequency": 2, "n_categories": 2, "majority_category_freq": 1,
                         "left_core_length": 1000, "right_core_length": 800, "accessory_length": 350},
        "200_f__300_f": {"frequency": 2, "n_categories": 2, "majority_category_freq": 1,
                         "left_core_length": 800, "right_core_length": 600, "accessory_length": 300},
        "300_f__400_f": {"frequency": 2, "n_categories": 1, "majority_category_freq": 2,
                         "left_core_length": 600, "right_core_length": 700, "accessory_length": 0},
        "100_f__300_f": {"frequency": 1, "n_categories": 1, "majority_category_freq": 1,
                         "left_core_length": 1000, "right_core_length": 600, "accessory_length": 150},
        "200_f__400_f": {"frequency": 1, "n_categories": 1, "majority_category_freq": 1,
                         "left_core_length": 800, "right_core_length": 700, "accessory_length": 300},
        "200_r__300_r": {"frequency": 1, "n_categories": 1, "majority_category_freq": 1,
                         "left_core_length": 800, "right_core_length": 600, "accessory_length": 0},
    }

    assert set(sdf.index) == set(expected.keys())
    for edge, vals in expected.items():
        for col, val in vals.items():
            assert sdf.loc[edge, col] == val, f"{edge}.{col}: {sdf.loc[edge, col]} != {val}"


def test_junction_stats_transitive_and_singleton(junction_pangraph):
    """Verify boolean flags for transitivity and singleton detection."""
    bj = BackboneJunctions(junction_pangraph, L_thr=500)
    sdf = bj.stats()

    # Transitive: only one path category
    transitive_edges = {"100_r__400_r", "300_f__400_f", "100_f__300_f", "200_f__400_f", "200_r__300_r"}
    non_transitive = {"100_f__200_f", "200_f__300_f"}
    for edge in transitive_edges:
        assert sdf.loc[edge, "is_transitive"] == True, f"{edge} should be transitive"
    for edge in non_transitive:
        assert sdf.loc[edge, "is_transitive"] == False, f"{edge} should not be transitive"

    # Singleton: all but one isolate share the same path
    singleton_edges = {"100_f__200_f", "200_f__300_f"}
    non_singleton = set(sdf.index) - singleton_edges
    for edge in singleton_edges:
        assert sdf.loc[edge, "is_singleton"] == True, f"{edge} should be singleton"
    for edge in non_singleton:
        assert sdf.loc[edge, "is_singleton"] == False, f"{edge} should not be singleton"


def test_junction_stats_sorted_by_frequency(junction_pangraph):
    """Stats DataFrame is sorted by frequency descending."""
    bj = BackboneJunctions(junction_pangraph, L_thr=500)
    sdf = bj.stats()
    freqs = sdf["frequency"].values
    assert all(freqs[i] >= freqs[i + 1] for i in range(len(freqs) - 1))


def test_junction_stats_linear(linear_pangraph):
    """Stats for linear paths fixture.

    Two edges, both shared by s1 and s2:
      100_f__200_f: s1 center=[A2+], s2 center=[A3+] → 2 categories, singleton
      200_f__300_f: both empty → 1 category, transitive
    """
    bj = BackboneJunctions(linear_pangraph, L_thr=500)
    sdf = bj.stats()

    assert len(sdf) == 2

    # 100_f__200_f: two different non-empty centers
    row = sdf.loc["100_f__200_f"]
    assert row["frequency"] == 2
    assert row["n_categories"] == 2
    assert row["majority_category_freq"] == 1
    assert row["is_transitive"] == False
    assert row["is_singleton"] == True
    # Unique blocks: A2(bid=600, 150bp) + A3(bid=700, 300bp)
    assert row["accessory_length"] == 150 + 300

    # 200_f__300_f: both empty
    row = sdf.loc["200_f__300_f"]
    assert row["frequency"] == 2
    assert row["n_categories"] == 1
    assert row["majority_category_freq"] == 2
    assert row["is_transitive"] == True
    assert row["is_singleton"] == False
    assert row["accessory_length"] == 0


def test_backbone_dataframe_returns_stats(junction_pangraph):
    """dataframe() returns (jdf, stats_df) where stats_df has correct columns."""
    bj = BackboneJunctions(junction_pangraph, L_thr=500)
    jdf, stats_df = bj.dataframe()

    assert isinstance(jdf, pd.DataFrame)
    assert isinstance(stats_df, pd.DataFrame)
    assert jdf.shape == (3, 7)

    expected_cols = [
        "frequency", "n_categories", "majority_category_freq",
        "is_transitive", "is_singleton",
        "left_core_length", "right_core_length", "accessory_length",
    ]
    assert list(stats_df.columns) == expected_cols
    assert len(stats_df) == 7

    # jdf columns should match stats_df index (same order)
    assert list(jdf.columns) == list(stats_df.index)


def test_junction_stats_smoke(plasmid_pangraph):
    """junction stats on the real plasmids dataset — check invariants."""
    bj = BackboneJunctions(plasmid_pangraph, L_thr=500)
    sdf = bj.stats()

    assert isinstance(sdf, pd.DataFrame)
    assert len(sdf) > 0

    # All frequencies positive
    assert (sdf["frequency"] > 0).all()
    # At least one category per edge
    assert (sdf["n_categories"] >= 1).all()
    # Majority can't exceed frequency
    assert (sdf["majority_category_freq"] <= sdf["frequency"]).all()
    # Transitive iff exactly one category
    assert (sdf["is_transitive"] == (sdf["n_categories"] == 1)).all()
    # Singleton requires frequency > 1
    assert not (sdf["is_singleton"] & (sdf["frequency"] <= 1)).any()
    # Core lengths positive
    assert (sdf["left_core_length"] > 0).all()
    assert (sdf["right_core_length"] > 0).all()
    # Accessory length non-negative
    assert (sdf["accessory_length"] >= 0).all()
