"""Shared test fixtures for pypangraph tests."""

import pytest
import pypangraph as pp


def _make_node(node_id, block_id, path_id, strand, start, end):
    return {
        "id": node_id,
        "block_id": block_id,
        "path_id": path_id,
        "strand": "+" if strand else "-",
        "position": [start, end],
    }


def _make_block(block_id, length, node_ids):
    return {
        "id": block_id,
        "consensus": "A" * length,
        "alignments": {
            str(nid): {"subs": [], "dels": [], "inss": []} for nid in node_ids
        },
    }


def _make_path(path_id, name, node_ids, tot_len, circular=True):
    return {
        "id": path_id,
        "nodes": node_ids,
        "tot_len": tot_len,
        "circular": circular,
        "name": name,
        "desc": None,
    }


def build_junction_pangraph_json():
    """Build a synthetic pangraph JSON with core and accessory blocks.

    Graph topology (all paths circular):
      s1: C1+ A1+ A2+ C2+ C3+ C4+
      s2: C1+ A1+ C2+ A3+ C3+ C4+
      s3: C1+ A2- C3+ C2+ A3+ C4+   (C2/C3 swapped = rearrangement)

    Core blocks (>=500bp, present once per strain):
      C1=100 (1000bp), C2=200 (800bp), C3=300 (600bp), C4=400 (700bp)
    Accessory blocks (<500bp):
      A1=500 (200bp, s1+s2), A2=600 (150bp, s1+s3), A3=700 (300bp, s2+s3)
    """
    # --- nodes ---
    # s1 (path 0): C1+(1) A1+(2) A2+(3) C2+(4) C3+(5) C4+(6)
    # s2 (path 1): C1+(7) A1+(8) C2+(9) A3+(10) C3+(11) C4+(12)
    # s3 (path 2): C1+(13) A2-(14) C3+(15) C2+(16) A3+(17) C4+(18)
    nodes = {
        # s1
        "1": _make_node(1, 100, 0, True, 0, 1000),
        "2": _make_node(2, 500, 0, True, 1000, 1200),
        "3": _make_node(3, 600, 0, True, 1200, 1350),
        "4": _make_node(4, 200, 0, True, 1350, 2150),
        "5": _make_node(5, 300, 0, True, 2150, 2750),
        "6": _make_node(6, 400, 0, True, 2750, 3450),
        # s2
        "7": _make_node(7, 100, 1, True, 0, 1000),
        "8": _make_node(8, 500, 1, True, 1000, 1200),
        "9": _make_node(9, 200, 1, True, 1200, 2000),
        "10": _make_node(10, 700, 1, True, 2000, 2300),
        "11": _make_node(11, 300, 1, True, 2300, 2900),
        "12": _make_node(12, 400, 1, True, 2900, 3600),
        # s3
        "13": _make_node(13, 100, 2, True, 0, 1000),
        "14": _make_node(14, 600, 2, False, 1000, 1150),
        "15": _make_node(15, 300, 2, True, 1150, 1750),
        "16": _make_node(16, 200, 2, True, 1750, 2550),
        "17": _make_node(17, 700, 2, True, 2550, 2850),
        "18": _make_node(18, 400, 2, True, 2850, 3550),
    }

    blocks = {
        "100": _make_block(100, 1000, [1, 7, 13]),
        "200": _make_block(200, 800, [4, 9, 16]),
        "300": _make_block(300, 600, [5, 11, 15]),
        "400": _make_block(400, 700, [6, 12, 18]),
        "500": _make_block(500, 200, [2, 8]),
        "600": _make_block(600, 150, [3, 14]),
        "700": _make_block(700, 300, [10, 17]),
    }

    paths = {
        "0": _make_path(0, "s1", [1, 2, 3, 4, 5, 6], 3450),
        "1": _make_path(1, "s2", [7, 8, 9, 10, 11, 12], 3600),
        "2": _make_path(2, "s3", [13, 14, 15, 16, 17, 18], 3550),
    }

    return {"paths": paths, "blocks": blocks, "nodes": nodes}


def build_linear_pangraph_json():
    """Build a synthetic pangraph JSON with linear (non-circular) paths.

    Graph topology (all paths linear):
      s1: A1+ C1+ A2+ C2+ C3+
      s2: C1+ A3+ C2+ C3+ A1+

    Core blocks (>=500bp, present once per strain):
      C1=100 (1000bp), C2=200 (800bp), C3=300 (600bp)
    Accessory blocks (<500bp):
      A1=500 (200bp, s1+s2), A2=600 (150bp, s1 only), A3=700 (300bp, s2 only)

    s1 has a leading accessory block (A1 before C1) and an internal junction.
    s2 has a trailing accessory block (A1 after C3) and an internal junction.
    """
    # s1 (path 0): A1+(1) C1+(2) A2+(3) C2+(4) C3+(5)
    # s2 (path 1): C1+(6) A3+(7) C2+(8) C3+(9) A1+(10)
    nodes = {
        # s1
        "1": _make_node(1, 500, 0, True, 0, 200),
        "2": _make_node(2, 100, 0, True, 200, 1200),
        "3": _make_node(3, 600, 0, True, 1200, 1350),
        "4": _make_node(4, 200, 0, True, 1350, 2150),
        "5": _make_node(5, 300, 0, True, 2150, 2750),
        # s2
        "6": _make_node(6, 100, 1, True, 0, 1000),
        "7": _make_node(7, 700, 1, True, 1000, 1300),
        "8": _make_node(8, 200, 1, True, 1300, 2100),
        "9": _make_node(9, 300, 1, True, 2100, 2700),
        "10": _make_node(10, 500, 1, True, 2700, 2900),
    }

    blocks = {
        "100": _make_block(100, 1000, [2, 6]),
        "200": _make_block(200, 800, [4, 8]),
        "300": _make_block(300, 600, [5, 9]),
        "500": _make_block(500, 200, [1, 10]),
        "600": _make_block(600, 150, [3]),
        "700": _make_block(700, 300, [7]),
    }

    paths = {
        "0": _make_path(0, "s1", [1, 2, 3, 4, 5], 2750, circular=False),
        "1": _make_path(1, "s2", [6, 7, 8, 9, 10], 2900, circular=False),
    }

    return {"paths": paths, "blocks": blocks, "nodes": nodes}


def _make_block_with_edits(block_id, consensus, node_edits):
    """Create a block with a custom consensus and per-node edits.

    Args:
        block_id: Block identifier.
        consensus: The consensus sequence string.
        node_edits: dict of node_id -> {"subs": [...], "dels": [...], "inss": [...]}
            or node_id -> None for identity (no edits).
    """
    alignments = {}
    for nid, edits in node_edits.items():
        if edits is None:
            edits = {"subs": [], "dels": [], "inss": []}
        alignments[str(nid)] = edits
    return {
        "id": block_id,
        "consensus": consensus,
        "alignments": alignments,
    }


def build_sequence_pangraph_json():
    """Minimal pangraph for testing junction sequence extraction.

    Two strains (circular), sharing the same junction in opposite orientations:
      s1: C1+(n1) A1+(n2) C2+(n3)    ← forward orientation
      s2: C2-(n4) A1-(n5) C1-(n6)    ← inverted orientation

    Block consensuses use distinct sequences for verification:
      C1 (bid=10, 6bp): "AAACCC" - backbone (len >= L_thr=4)
      C2 (bid=20, 6bp): "GGGAAA" - backbone
      A1 (bid=30, 3bp): "TTC"    - accessory (len < L_thr=4)

    Node 6 (C1 in s2) has a substitution: pos=0 A->T, giving "TAACCC".
    This makes the two isolates' sequences distinguishable.
    """
    nodes = {
        # s1
        "1": _make_node(1, 10, 0, True, 0, 6),
        "2": _make_node(2, 30, 0, True, 6, 9),
        "3": _make_node(3, 20, 0, True, 9, 15),
        # s2
        "4": _make_node(4, 20, 1, False, 0, 6),
        "5": _make_node(5, 30, 1, False, 6, 9),
        "6": _make_node(6, 10, 1, False, 9, 15),
    }

    blocks = {
        "10": _make_block_with_edits(10, "AAACCC", {
            1: None,
            6: {"subs": [{"pos": 0, "alt": "T"}], "dels": [], "inss": []},
        }),
        "20": _make_block_with_edits(20, "GGGAAA", {
            3: None,
            4: None,
        }),
        "30": _make_block_with_edits(30, "TTC", {
            2: None,
            5: None,
        }),
    }

    paths = {
        "0": _make_path(0, "s1", [1, 2, 3], 15),
        "1": _make_path(1, "s2", [4, 5, 6], 15),
    }

    return {"paths": paths, "blocks": blocks, "nodes": nodes}


@pytest.fixture
def sequence_pangraph():
    """A minimal Pangraph for testing junction sequence extraction."""
    return pp.Pangraph(build_sequence_pangraph_json())


@pytest.fixture
def junction_pangraph():
    """A synthetic Pangraph with core/accessory blocks and non-trivial junctions."""
    return pp.Pangraph(build_junction_pangraph_json())


@pytest.fixture
def linear_pangraph():
    """A synthetic Pangraph with linear (non-circular) paths."""
    return pp.Pangraph(build_linear_pangraph_json())
