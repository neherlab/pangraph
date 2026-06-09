"""Tests for the GFA writer and the streamlined junction-aware GFA export."""

import math
import re

import pytest

from pypangraph.export import GFA, junction_context_gfa
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


def test_gfa_write_minimal(tmp_path):
    """GFA.write emits H/S/L lines with LN/DP tags and correct orientations."""
    segments = {"a": 100, "b": 50}
    links = {("a", True, "b", False)}
    depths = {"a": 3}
    out = tmp_path / "tiny.gfa"
    GFA(segments, links, depths).write(out)

    lines = out.read_text().splitlines()
    assert lines[0] == "H\tVN:Z:1.0"
    s_lines = {ln.split("\t")[1]: ln for ln in lines if ln.startswith("S")}
    assert "LN:i:100" in s_lines["a"] and "DP:f:3" in s_lines["a"]
    assert "LN:i:50" in s_lines["b"] and "DP:f:" not in s_lines["b"]
    (link,) = [ln for ln in lines if ln.startswith("L")]
    assert link == "L\ta\t+\tb\t-\t0M"


def _segment_lengths(gfa_path):
    """Parse a GFA file into a dict of segment name -> emitted LN:i: value."""
    lengths = {}
    for ln in gfa_path.read_text().splitlines():
        if not ln.startswith("S"):
            continue
        fields = ln.split("\t")
        (tag,) = [f for f in fields if f.startswith("LN:i:")]
        lengths[fields[1]] = int(tag.removeprefix("LN:i:"))
    return lengths


def test_gfa_write_length_transform(tmp_path):
    """length_transform rescales LN:i:, rounds to nearest int, clamps to >= 1."""
    gfa = GFA({"a": 100, "b": 1000, "c": 5}, set(), {})

    # Linear rescale: round(l / 10) -> exact for 100 and 1000; 5 -> round(0.5)=0
    # clamps up to the 1 floor.
    out = tmp_path / "scaled.gfa"
    gfa.write(out, length_transform=lambda length: length / 10)
    assert _segment_lengths(out) == {"a": 10, "b": 100, "c": 1}

    # log-based scaling also rounds to the nearest integer.
    out_log = tmp_path / "log.gfa"
    gfa.write(out_log, length_transform=lambda length: 100 * math.log(length))
    assert _segment_lengths(out_log) == {
        name: max(1, round(100 * math.log(length)))
        for name, length in gfa.segments.items()
    }

    # Default path is unaffected: real lengths.
    out_real = tmp_path / "real.gfa"
    gfa.write(out_real)
    assert _segment_lengths(out_real) == {"a": 100, "b": 1000, "c": 5}


def test_consensus_gfa_structure(junction_pangraph):
    """Consensus export: core anchors single & un-prefixed, accessory prefixed,
    links well-formed, and depths reflect isolate coverage."""
    bj = BackboneJunctions(junction_pangraph, L_thr=500)
    gfa, prefix_map = junction_context_gfa(bj, scaffold="consensus")
    segments, links, depths = gfa.segments, gfa.links, gfa.depths

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


def test_all_scaffold_is_superset_of_consensus(junction_pangraph):
    """scaffold='all' keeps every junction, a superset of the consensus scaffold."""
    bj = BackboneJunctions(junction_pangraph, L_thr=500)
    gfa_c, pmap_c = junction_context_gfa(bj, scaffold="consensus")
    gfa_a, pmap_a = junction_context_gfa(bj, scaffold="all")

    edges_c = set(pmap_c.values())
    edges_a = set(pmap_a.values())
    assert edges_c < edges_a  # strict superset (s3 rearrangement edges added)
    assert len(gfa_a.links) >= len(gfa_c.links)


def test_consensus_scaffold_follows_dominant_synteny(junction_pangraph):
    """Per-edge majority keeps the s1/s2 backbone, dropping the s3 rearrangement."""
    bj = BackboneJunctions(junction_pangraph, L_thr=500)
    _, prefix_map = junction_context_gfa(bj, scaffold="consensus")
    assert set(prefix_map.values()) == DOMINANT_EDGES


def test_reference_scaffold_uses_that_genomes_edges(junction_pangraph):
    """scaffold=<isolate> keeps that genome's own core edges (s3 is rearranged)."""
    bj = BackboneJunctions(junction_pangraph, L_thr=500)
    _, prefix_map = junction_context_gfa(bj, scaffold="s3")
    edges = set(prefix_map.values())
    # s3 visits C1-C3-C2-C4, so its core edges differ from the consensus backbone.
    assert edges != DOMINANT_EDGES
    assert "100_f__300_f" in edges


def test_unknown_scaffold_raises(junction_pangraph):
    """A scaffold that is neither 'consensus'/'all' nor a known isolate raises."""
    bj = BackboneJunctions(junction_pangraph, L_thr=500)
    with pytest.raises(ValueError, match="unknown scaffold isolate"):
        junction_context_gfa(bj, scaffold="not_a_genome")


def test_junction_context_gfa_runs_on_real_graph(plasmid_pangraph):
    """Smoke test on a real plasmid graph: non-empty, all links valid."""
    bj = BackboneJunctions(plasmid_pangraph, L_thr=500)
    gfa, prefix_map = junction_context_gfa(bj, scaffold="consensus")
    assert gfa.segments and gfa.links and prefix_map
    for a, _, b, _ in gfa.links:
        assert a in gfa.segments and b in gfa.segments
