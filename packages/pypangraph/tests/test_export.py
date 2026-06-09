"""Tests for the GFA writer and the streamlined junction-aware GFA export."""

import re

import pytest

from pypangraph.export import streamlined_gfa, write_gfa
from pypangraph.junctions import BackboneJunctions

ACCESSORY_RE = re.compile(r"^J\d+__\d+$")

# The dominant backbone shared by s1 and s2 (C1-C2-C3-C4, circular); s3 carries
# the rearranged C1-C3-C2-C4 instead. See conftest.build_junction_pangraph_json.
DOMINANT_EDGES = {
    "100_f__200_f",
    "200_f__300_f",
    "300_f__400_f",
    "100_r__400_r",
}


def test_write_gfa_minimal(tmp_path):
    """write_gfa emits H/S/L lines with LN/DP tags and correct orientations."""
    segments = {"a": 100, "b": 50}
    links = {("a", True, "b", False)}
    depths = {"a": 3}
    out = tmp_path / "tiny.gfa"
    write_gfa(out, segments, links, depths)

    lines = out.read_text().splitlines()
    assert lines[0] == "H\tVN:Z:1.0"
    s_lines = {ln.split("\t")[1]: ln for ln in lines if ln.startswith("S")}
    assert "LN:i:100" in s_lines["a"] and "DP:f:3" in s_lines["a"]
    assert "LN:i:50" in s_lines["b"] and "DP:f:" not in s_lines["b"]
    (link,) = [ln for ln in lines if ln.startswith("L")]
    assert link == "L\ta\t+\tb\t-\t0M"


def test_streamlined_consensus_structure(junction_pangraph):
    """Consensus export: core anchors single & un-prefixed, accessory prefixed,
    links well-formed, and depths reflect isolate coverage."""
    bj = BackboneJunctions(junction_pangraph, L_thr=500)
    segments, links, depths, prefix_map = streamlined_gfa(bj, order="consensus")

    core_ids = {"100", "200", "300", "400"}
    # Every core block appears once, un-prefixed.
    assert core_ids <= set(segments)
    # Non-core segments are all per-junction accessory copies.
    for name in segments:
        if name not in core_ids:
            assert ACCESSORY_RE.match(name), name

    # Links only reference existing segments.
    for a, _, b, _ in links:
        assert a in segments and b in segments

    # Core anchors are present in every isolate, so all share the same depth (3)
    # regardless of which rearranged edges the consensus scaffold drops.
    for cid in core_ids:
        assert depths[cid] == 3
    # A1 (id 500) sits in the C1-C2 junction of s1 and s2 -> depth 2.
    a1 = next(n for n in segments if n.endswith("__500"))
    assert depths[a1] == 2

    # Every kept junction maps to a real edge of the graph.
    for edge_str in prefix_map.values():
        assert edge_str in bj


def test_streamlined_all_is_superset_of_consensus(junction_pangraph):
    """order='all' keeps every junction, a superset of the consensus scaffold."""
    bj = BackboneJunctions(junction_pangraph, L_thr=500)
    _, links_c, _, pmap_c = streamlined_gfa(bj, order="consensus")
    _, links_a, _, pmap_a = streamlined_gfa(bj, order="all")

    edges_c = set(pmap_c.values())
    edges_a = set(pmap_a.values())
    assert edges_c < edges_a  # strict superset (s3 rearrangement edges added)
    assert len(links_a) >= len(links_c)


def test_consensus_scaffold_follows_dominant_synteny(junction_pangraph):
    """Per-edge majority keeps the s1/s2 backbone, dropping the s3 rearrangement."""
    bj = BackboneJunctions(junction_pangraph, L_thr=500)
    _, _, _, prefix_map = streamlined_gfa(bj, order="consensus")
    assert set(prefix_map.values()) == DOMINANT_EDGES


def test_reference_scaffold_uses_that_genomes_edges(junction_pangraph):
    """order=<isolate> keeps that genome's own core edges (s3 is rearranged)."""
    bj = BackboneJunctions(junction_pangraph, L_thr=500)
    _, _, _, prefix_map = streamlined_gfa(bj, order="s3")
    edges = set(prefix_map.values())
    # s3 visits C1-C3-C2-C4, so its core edges differ from the consensus backbone.
    assert edges != DOMINANT_EDGES
    assert "100_f__300_f" in edges


def test_unknown_order_raises(junction_pangraph):
    """An order that is neither 'consensus'/'all' nor a known isolate raises."""
    bj = BackboneJunctions(junction_pangraph, L_thr=500)
    with pytest.raises(ValueError, match="unknown reference isolate"):
        streamlined_gfa(bj, order="not_a_genome")


def test_streamlined_runs_on_real_graph(plasmid_pangraph):
    """Smoke test on a real plasmid graph: non-empty, all links valid."""
    bj = BackboneJunctions(plasmid_pangraph, L_thr=500)
    segments, links, depths, prefix_map = streamlined_gfa(bj, order="consensus")
    assert segments and links and prefix_map
    for a, _, b, _ in links:
        assert a in segments and b in segments
