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
      s1: c5+ C1+ A2+ C2+ C3+
      s2: C1+ A3+ C2+ C3+ c5+

    Backbone blocks (core AND >=500bp, present once per strain):
      C1=100 (1000bp), C2=200 (800bp), C3=300 (600bp)
    Accessory blocks (present in only one strain):
      A2=600 (150bp, s1 only), A3=700 (300bp, s2 only)
    Core-but-not-backbone block (lowercase 'c' to flag it):
      c5=500 (200bp, present once in both strains)

    c5 is strictly a `core` block (it occurs exactly once in every strain), but its
    200bp length is below the backbone threshold used by the junction/MSU analyses
    (core AND len>=L_thr), so it cannot be a backbone block and is handled as accessory
    content. It is the same block at a leading position in s1 and a trailing position in
    s2, exercising terminal junctions at both ends.
    """
    # s1 (path 0): c5+(1) C1+(2) A2+(3) C2+(4) C3+(5)
    # s2 (path 1): C1+(6) A3+(7) C2+(8) C3+(9) c5+(10)
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
      s2: C2-(n4) A2-(n5) C1-(n6)    ← inverted orientation

    Block consensuses use distinct sequences for verification:
      C1 (bid=10, 6bp): "AAACCC" - backbone (len >= L_thr=4)
      C2 (bid=20, 6bp): "GGGAAA" - backbone
      A1 (bid=30, 3bp): "TTC"    - accessory (s1 only, len < L_thr=4)
      A2 (bid=40, 3bp): "ACG"    - accessory (s2 only, len < L_thr=4)

    A1 and A2 are distinct strain-private accessory blocks with different consensus
    sequences, so neither junction center block is `core`. Node 6 (C1 in s2) has a
    substitution: pos=0 A->T, giving "TAACCC", which makes the isolates' flanking
    sequences distinguishable too.
    """
    nodes = {
        # s1
        "1": _make_node(1, 10, 0, True, 0, 6),
        "2": _make_node(2, 30, 0, True, 6, 9),
        "3": _make_node(3, 20, 0, True, 9, 15),
        # s2
        "4": _make_node(4, 20, 1, False, 0, 6),
        "5": _make_node(5, 40, 1, False, 6, 9),
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
        }),
        "40": _make_block_with_edits(40, "ACG", {
            5: None,
        }),
    }

    paths = {
        "0": _make_path(0, "s1", [1, 2, 3], 15),
        "1": _make_path(1, "s2", [4, 5, 6], 15),
    }

    return {"paths": paths, "blocks": blocks, "nodes": nodes}


# Consensus sequences for build_inversion_pangraph_json, keyed by block-id string.
# Distinct and non-palindromic so junction sequence extraction (incl. reverse-complement)
# can be verified. Lengths are variable and <20bp; the seven core blocks (10-70) are >10bp.
INVERSION_CONS = {
    "10": "ACGTTGCAACCA",      # C1, 12bp
    "20": "TTGGAACCGGTTAC",    # C2, 14bp
    "30": "GATTACAGGCT",       # C3, 11bp
    "40": "CCAGTACGTGACATCA",  # C4, 16bp
    "50": "ACACGTGTACGTA",     # C5, 13bp
    "60": "TGTCATGCAATGCAT",   # C6, 15bp
    "70": "GGATCCGAATTCAGTCA",  # C7, 17bp
    "80": "ACGTGA",            # A1, 6bp (accessory)
    "90": "TTCAGGCA",          # A2, 8bp (accessory)
}


def build_inversion_pangraph_json():
    """Circular pangraph exercising reverse-complement, a single inversion, and mergers.

    Three strains (all circular). C5/C6/C7 sit next to and co-oriented with C1/C2/C3, so
    each pair always travels together and forms a merger (MSU):
      s1 (reference):       C1+ C5+ A1+ C2+ C6+ C3+ C7+ A2- C4+
      s2 (whole-genome RC): C4- A2+ C7- C3- C6- C2- A1- C5- C1-
      s3 (single inversion): C1+ C5+ C2+ C6+ C7- C3- C4+   (no accessory)

    s2 is the exact reverse complement of s1 (order reversed, every strand flipped), so it
    carries the same backbone edges as s1 (Edge equality is orientation-canonical) and must
    co-orient back onto s1. s3 inverts the contiguous `C3 C7` segment: the internal C3-C7
    adjacency is preserved (so {C3,C7} still merge) while C6-C3 and C7-C4 are broken, leaving
    {C3,C7} as a separate, invertible MSU. A2 is on the `-` strand in s1 so the junction
    sequence extraction exercises the center reverse-complement path. See INVERSION_CONS for
    block consensuses (variable <20bp lengths; core blocks >10bp; tested at L_thr=10).

    Core blocks (once per strain): C1=10, C2=20, C3=30, C4=40, C5=50, C6=60, C7=70.
    Accessory blocks (s1+s2 only): A1=80, A2=90.
    """
    nodes = {
        # s1: C1+ C5+ A1+ C2+ C6+ C3+ C7+ A2- C4+
        "1": _make_node(1, 10, 0, True, 0, 12),
        "2": _make_node(2, 50, 0, True, 12, 25),
        "3": _make_node(3, 80, 0, True, 25, 31),
        "4": _make_node(4, 20, 0, True, 31, 45),
        "5": _make_node(5, 60, 0, True, 45, 60),
        "6": _make_node(6, 30, 0, True, 60, 71),
        "7": _make_node(7, 70, 0, True, 71, 88),
        "8": _make_node(8, 90, 0, False, 88, 96),
        "9": _make_node(9, 40, 0, True, 96, 112),
        # s2: C4- A2+ C7- C3- C6- C2- A1- C5- C1-
        "10": _make_node(10, 40, 1, False, 0, 16),
        "11": _make_node(11, 90, 1, True, 16, 24),
        "12": _make_node(12, 70, 1, False, 24, 41),
        "13": _make_node(13, 30, 1, False, 41, 52),
        "14": _make_node(14, 60, 1, False, 52, 67),
        "15": _make_node(15, 20, 1, False, 67, 81),
        "16": _make_node(16, 80, 1, False, 81, 87),
        "17": _make_node(17, 50, 1, False, 87, 100),
        "18": _make_node(18, 10, 1, False, 100, 112),
        # s3: C1+ C5+ C2+ C6+ C7- C3- C4+
        "19": _make_node(19, 10, 2, True, 0, 12),
        "20": _make_node(20, 50, 2, True, 12, 25),
        "21": _make_node(21, 20, 2, True, 25, 39),
        "22": _make_node(22, 60, 2, True, 39, 54),
        "23": _make_node(23, 70, 2, False, 54, 71),
        "24": _make_node(24, 30, 2, False, 71, 82),
        "25": _make_node(25, 40, 2, True, 82, 98),
    }

    block_nodes = {
        "10": [1, 18, 19],
        "20": [4, 15, 21],
        "30": [6, 13, 24],
        "40": [9, 10, 25],
        "50": [2, 17, 20],
        "60": [5, 14, 22],
        "70": [7, 12, 23],
        "80": [3, 16],
        "90": [8, 11],
    }
    blocks = {
        bid: _make_block_with_edits(
            int(bid), INVERSION_CONS[bid], {nid: None for nid in nids}
        )
        for bid, nids in block_nodes.items()
    }

    paths = {
        "0": _make_path(0, "s1", [1, 2, 3, 4, 5, 6, 7, 8, 9], 112),
        "1": _make_path(1, "s2", [10, 11, 12, 13, 14, 15, 16, 17, 18], 112),
        "2": _make_path(2, "s3", [19, 20, 21, 22, 23, 24, 25], 98),
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


@pytest.fixture
def inversion_pangraph():
    """A circular Pangraph with a reverse-complemented strain, a single inversion, and mergers."""
    return pp.Pangraph(build_inversion_pangraph_json())
