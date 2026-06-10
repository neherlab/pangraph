"""Tests for the BackboneJunctions analyses: positions(), sequences(), stats()."""

import pandas as pd

from pypangraph.junctions import BackboneJunctions


# --- positions tests ---


def test_junction_positions_forward_strand(junction_pangraph):
    """Verify positions for forward-strand junctions in s1.

    s1 path: C1+(0,1000) A1+(1000,1200) A2+(1200,1350) C2+(1350,2150) C3+(2150,2750) C4+(2750,3450)
    """
    bj = BackboneJunctions(junction_pangraph, L_thr=500)
    pos = bj.positions()

    # Edge C1+→C2+ in s1: both flanks on forward strand → strand=True
    row = pos.loc[("100_f__200_f", "s1")]
    assert row["strand"]
    assert row["left_start"] == 0
    assert row["left_end"] == 1000
    assert row["right_start"] == 1350
    assert row["right_end"] == 2150

    # Edge C2+→C3+ in s1: empty junction, strand=True
    row = pos.loc[("200_f__300_f", "s1")]
    assert row["strand"]
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
    bj = BackboneJunctions(junction_pangraph, L_thr=500)
    pos = bj.positions()

    # In s1: junction is inverted, left=C4, right=C1
    row = pos.loc[("100_r__400_r", "s1")]
    assert not row["strand"]
    assert row["left_start"] == 2750  # C4 start
    assert row["left_end"] == 3450  # C4 end
    assert row["right_start"] == 0  # C1 start
    assert row["right_end"] == 1000  # C1 end


def test_junction_positions_rearranged_strain(junction_pangraph):
    """Verify positions for s3 which has a rearrangement (C2/C3 swapped).

    s3 path: C1+(0,1000) A2-(1000,1150) C3+(1150,1750) C2+(1750,2550) A3+(2550,2850) C4+(2850,3550)
    """
    bj = BackboneJunctions(junction_pangraph, L_thr=500)
    pos = bj.positions()

    # Edge C1+→C3+: forward in s3
    row = pos.loc[("100_f__300_f", "s3")]
    assert row["strand"]
    assert row["left_start"] == 0  # C1
    assert row["left_end"] == 1000
    assert row["right_start"] == 1150  # C3
    assert row["right_end"] == 1750

    # Edge "200_r__300_r" = Edge(C2-, C3-): in s3 both are on + strand → inverted
    # Swapped: left=C3(1150,1750), right=C2(1750,2550)
    row = pos.loc[("200_r__300_r", "s3")]
    assert not row["strand"]
    assert row["left_start"] == 1150  # C3
    assert row["left_end"] == 1750
    assert row["right_start"] == 1750  # C2
    assert row["right_end"] == 2550

    # Edge C2+→C4+: forward in s3, A3 between them
    row = pos.loc[("200_f__400_f", "s3")]
    assert row["strand"]
    assert row["left_start"] == 1750  # C2
    assert row["left_end"] == 2550
    assert row["right_start"] == 2850  # C4
    assert row["right_end"] == 3550


def test_junction_positions_shape(junction_pangraph):
    """Position DataFrame has one row per (isolate, edge) with non-NaN junction length."""
    bj = BackboneJunctions(junction_pangraph, L_thr=500)
    pos = bj.positions()

    # one position row per (iso, edge) carrying that junction
    n_junctions = bj.stats()["n_isolates"].sum()
    assert len(pos) == n_junctions
    assert list(pos.columns) == [
        "left_start",
        "left_end",
        "right_start",
        "right_end",
        "strand",
    ]


def test_junction_positions_linear(linear_pangraph):
    """Positions work correctly for linear paths.

    s1: c5+(0,200) C1+(200,1200) A2+(1200,1350) C2+(1350,2150) C3+(2150,2750)
    s2: C1+(0,1000) A3+(1000,1300) C2+(1300,2100) C3+(2100,2700) c5+(2700,2900)
    """
    bj = BackboneJunctions(linear_pangraph, L_thr=500)
    pos = bj.positions()

    # terminal junctions (leading c5 in s1, trailing c5 in s2) have no flanking
    # edge, so only the two edge-bearing junctions appear
    assert set(pos.index.get_level_values("edge")) == {"100_f__200_f", "200_f__300_f"}

    # s1: C1→C2 has A2 in between
    row = pos.loc[("100_f__200_f", "s1")]
    assert row["left_start"] == 200  # C1
    assert row["left_end"] == 1200
    assert row["right_start"] == 1350  # C2
    assert row["right_end"] == 2150

    # s2: C1→C2 has A3 in between
    row = pos.loc[("100_f__200_f", "s2")]
    assert row["left_start"] == 0  # C1
    assert row["left_end"] == 1000
    assert row["right_start"] == 1300  # C2
    assert row["right_end"] == 2100


def test_junction_positions_smoke(plasmid_pangraph):
    """positions() runs on the real plasmids dataset without errors."""
    bj = BackboneJunctions(plasmid_pangraph, L_thr=500)
    pos = bj.positions()

    assert isinstance(pos, pd.DataFrame)
    n_junctions = bj.stats()["n_isolates"].sum()
    assert len(pos) == n_junctions
    assert list(pos.columns) == [
        "left_start",
        "left_end",
        "right_start",
        "right_end",
        "strand",
    ]


# --- sequence extraction tests ---


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

    s2: C2-(n4) A2-(n5) C1-(n6) — inverted for edge 10_f__20_f
    After inversion: C1+(n6) A2+(n5) C2+(n4)
    C1 n6 has sub pos=0 A->T, so seq = TAACCC
    A2 n5: ACG (no edits)
    C2 n4: GGGAAA (no edits)
    Expected: TAACCC + ACG + GGGAAA = TAACCCACGGGGAAA
    """
    bj = BackboneJunctions(sequence_pangraph, L_thr=4)
    records = bj.sequences("10_f__20_f")
    by_iso = {r.id: str(r.seq) for r in records}
    assert by_iso["s2"] == "TAACCCACGGGGAAA"


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
    stats_df = bj.stats()

    # Test on the most frequent edge
    top_edge = stats_df.index[0]
    top_freq = stats_df.loc[top_edge, "n_isolates"]
    records = bj.sequences(top_edge)
    assert len(records) == top_freq
    for r in records:
        assert len(r.seq) > 0
        assert r.description == top_edge


# --- junction stats tests ---


def test_junction_stats_values(junction_pangraph):
    """Verify all per-edge statistics against expected values.

    See plan for derivation of expected values.
    """
    bj = BackboneJunctions(junction_pangraph, L_thr=500)
    sdf = bj.stats()

    expected = {
        "100_r__400_r": {
            "n_isolates": 3,
            "n_non_empty": 0,
            "n_categories": 1,
            "n_majority_category": 3,
            "left_core_length": 1000,
            "right_core_length": 700,
            "accessory_length": 0,
        },
        "100_f__200_f": {
            "n_isolates": 2,
            "n_non_empty": 2,
            "n_categories": 2,
            "n_majority_category": 1,
            "left_core_length": 1000,
            "right_core_length": 800,
            "accessory_length": 350,
        },
        "200_f__300_f": {
            "n_isolates": 2,
            "n_non_empty": 1,
            "n_categories": 2,
            "n_majority_category": 1,
            "left_core_length": 800,
            "right_core_length": 600,
            "accessory_length": 300,
        },
        "300_f__400_f": {
            "n_isolates": 2,
            "n_non_empty": 0,
            "n_categories": 1,
            "n_majority_category": 2,
            "left_core_length": 600,
            "right_core_length": 700,
            "accessory_length": 0,
        },
        "100_f__300_f": {
            "n_isolates": 1,
            "n_non_empty": 1,
            "n_categories": 1,
            "n_majority_category": 1,
            "left_core_length": 1000,
            "right_core_length": 600,
            "accessory_length": 150,
        },
        "200_f__400_f": {
            "n_isolates": 1,
            "n_non_empty": 1,
            "n_categories": 1,
            "n_majority_category": 1,
            "left_core_length": 800,
            "right_core_length": 700,
            "accessory_length": 300,
        },
        "200_r__300_r": {
            "n_isolates": 1,
            "n_non_empty": 0,
            "n_categories": 1,
            "n_majority_category": 1,
            "left_core_length": 800,
            "right_core_length": 600,
            "accessory_length": 0,
        },
    }

    assert set(sdf.index) == set(expected.keys())
    for edge, vals in expected.items():
        for col, val in vals.items():
            assert sdf.loc[edge, col] == val, (
                f"{edge}.{col}: {sdf.loc[edge, col]} != {val}"
            )


def test_junction_stats_transitive_and_singleton(junction_pangraph):
    """Verify boolean flags for transitivity and singleton detection."""
    bj = BackboneJunctions(junction_pangraph, L_thr=500)
    sdf = bj.stats()

    # Transitive: only one path category
    transitive_edges = {
        "100_r__400_r",
        "300_f__400_f",
        "100_f__300_f",
        "200_f__400_f",
        "200_r__300_r",
    }
    non_transitive = {"100_f__200_f", "200_f__300_f"}
    for edge in transitive_edges:
        assert sdf.loc[edge, "is_transitive"], f"{edge} should be transitive"
    for edge in non_transitive:
        assert not sdf.loc[edge, "is_transitive"], f"{edge} should not be transitive"

    # Singleton: all but one isolate share the same path
    singleton_edges = {"100_f__200_f", "200_f__300_f"}
    non_singleton = set(sdf.index) - singleton_edges
    for edge in singleton_edges:
        assert sdf.loc[edge, "is_singleton"], f"{edge} should be singleton"
    for edge in non_singleton:
        assert not sdf.loc[edge, "is_singleton"], f"{edge} should not be singleton"


def test_junction_stats_sorted_by_n_isolates(junction_pangraph):
    """Stats DataFrame is sorted by `n_isolates` descending."""
    bj = BackboneJunctions(junction_pangraph, L_thr=500)
    sdf = bj.stats()
    freqs = sdf["n_isolates"].values
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
    # terminal junctions (no flanking edge) are excluded
    assert set(sdf.index) == {"100_f__200_f", "200_f__300_f"}

    # 100_f__200_f: two different non-empty centers
    row = sdf.loc["100_f__200_f"]
    assert row["n_isolates"] == 2
    assert row["n_non_empty"] == 2
    assert row["n_categories"] == 2
    assert row["n_majority_category"] == 1
    assert not row["is_transitive"]
    assert row["is_singleton"]
    # Unique blocks: A2(bid=600, 150bp) + A3(bid=700, 300bp)
    assert row["accessory_length"] == 150 + 300

    # 200_f__300_f: both empty
    row = sdf.loc["200_f__300_f"]
    assert row["n_isolates"] == 2
    assert row["n_non_empty"] == 0
    assert row["n_categories"] == 1
    assert row["n_majority_category"] == 2
    assert row["is_transitive"]
    assert not row["is_singleton"]
    assert row["accessory_length"] == 0


def test_backbone_stats_columns(junction_pangraph):
    """stats() returns a DataFrame with the documented column set."""
    bj = BackboneJunctions(junction_pangraph, L_thr=500)
    stats_df = bj.stats()

    assert isinstance(stats_df, pd.DataFrame)
    assert list(stats_df.columns) == [
        "n_isolates",
        "n_non_empty",
        "n_categories",
        "n_majority_category",
        "is_transitive",
        "is_singleton",
        "left_core_length",
        "right_core_length",
        "accessory_length",
    ]
    assert len(stats_df) == 7


def test_junction_stats_smoke(plasmid_pangraph):
    """junction stats on the real plasmids dataset — check invariants."""
    bj = BackboneJunctions(plasmid_pangraph, L_thr=500)
    sdf = bj.stats()

    assert isinstance(sdf, pd.DataFrame)
    assert len(sdf) > 0

    # All isolate counts positive
    assert (sdf["n_isolates"] > 0).all()
    # Non-empty count is bounded by `n_isolates`
    assert (sdf["n_non_empty"] >= 0).all()
    assert (sdf["n_non_empty"] <= sdf["n_isolates"]).all()
    # At least one category per edge
    assert (sdf["n_categories"] >= 1).all()
    # Majority can't exceed `n_isolates`
    assert (sdf["n_majority_category"] <= sdf["n_isolates"]).all()
    # Transitive iff exactly one category
    assert (sdf["is_transitive"] == (sdf["n_categories"] == 1)).all()
    # Singleton requires `n_isolates` > 1
    assert not (sdf["is_singleton"] & (sdf["n_isolates"] <= 1)).any()
    # Core lengths positive
    assert (sdf["left_core_length"] > 0).all()
    assert (sdf["right_core_length"] > 0).all()
    # Accessory length non-negative
    assert (sdf["accessory_length"] >= 0).all()
