"""Junction-context GFA decomposition of a pangenome graph.

Unlike the standard whole-graph GFA export, this disentangles the graph by
**junction context**: blocks are paralog-split and de-duplicated per junction
(the same block in two different junctions becomes two segments; shared core
anchors stay single), and the topology is reduced to a chosen **core synteny**
scaffold so the result is a clean, walkable graph instead of a tangle.

Core flank blocks are emitted once globally (segment name = block id); accessory
blocks get an enumerated per-junction prefix (``J{n}__{block_id}``) so blocks
shared across junctions are emitted once per junction. Each isolate's junction is
co-oriented to its canonical edge direction so the junction forms a clean bubble
between its two core anchors.
"""

from __future__ import annotations

from collections import Counter, defaultdict

from ..minimal_synteny_units import core_paths
from .gfa import GFA


def _iso_core_edges(pan, L_thr: int) -> dict[str, frozenset[str]]:
    """Map each isolate to the frozenset of canonical core-edge ids on its backbone.

    ``core_paths`` purifies each genome to core blocks >= ``L_thr``; the canonical
    ``Edge.to_str_id`` makes a reverse traversal yield identical ids. This works
    unchanged for multi-replicon genomes: edges from several disjoint core cycles
    just land in the same set.
    """
    return {
        iso: frozenset(e.to_str_id() for e in w.edges())
        for iso, w in core_paths(pan, L_thr).items()
    }


def _consensus_edge_set(iso_edges: dict[str, frozenset[str]]) -> set[str]:
    """Core edges present in a strict majority of isolate backbones.

    This is the single tunable consensus-policy point: change the predicate here
    (e.g. weighted or threshold-based voting) without touching the rest of the
    export.
    """
    n = len(iso_edges)
    counts = Counter(e for edges in iso_edges.values() for e in edges)
    return {e for e, c in counts.items() if c > n / 2}


def _scaffold_edges(bj, scaffold: str) -> list[str]:
    """Return the canonical edge ids to keep for a ``scaffold`` policy.

    The result is sorted for deterministic, stable ``J{n}`` numbering. ``"all"``
    keeps every junction (equivalently, the union of all isolate core-edge sets);
    ``"consensus"`` keeps the per-edge majority; otherwise ``scaffold`` is an
    isolate name and its own core-edge set is used. Edges absent from ``bj`` are
    dropped.
    """
    available = set(bj.edges())
    if scaffold == "all":
        return sorted(available)

    iso_edges = _iso_core_edges(bj.pan, bj.L_thr)
    if scaffold == "consensus":
        chosen = _consensus_edge_set(iso_edges)
    elif scaffold in iso_edges:
        chosen = iso_edges[scaffold]
    else:
        raise ValueError(
            f"unknown scaffold isolate {scaffold!r}; "
            f"expected 'consensus', 'all', or one of {sorted(iso_edges)}"
        )
    return sorted(chosen & available)


def junction_context_gfa(bj, scaffold: str = "consensus") -> tuple[GFA, dict[str, str]]:
    """Build a junction-context GFA decomposition of a graph.

    Args:
        bj: A ``BackboneJunctions`` object.
        scaffold: Which core-synteny scaffold defines the kept junctions.
            ``"consensus"`` (default) keeps the per-edge majority backbone; an
            isolate name uses that genome's own core edges; ``"all"`` keeps every
            junction (tangled escape hatch).

    Returns:
        A tuple ``(gfa, prefix_map)``:
        - gfa: a :class:`~pypangraph.export.GFA` whose segments (name -> length in
          bp), links and per-segment depths (number of distinct isolates
          traversing each segment across the kept junctions) describe the
          decomposition. Serialize it with ``gfa.write(path)``.
        - prefix_map: dict[str, str]  ``"J{n}"`` -> canonical edge string id;
          decomposition provenance the caller may persist (e.g. as a TSV).
    """
    bdf = bj.block_stats
    kept = _scaffold_edges(bj, scaffold)

    segments = {}
    links = set()
    prefix_map = {}
    seg_isolates = defaultdict(set)  # segment name -> set of isolates traversing it

    for n, edge_str in enumerate(kept):
        prefix = f"J{n}"
        prefix_map[prefix] = edge_str

        for iso, junction in bj[edge_str].items():
            jc = junction.to_canonical()
            core_ids = {jc.left.id, jc.right.id}

            named = []
            for ob in jc.oriented_blocks():
                name = ob.id if ob.id in core_ids else f"{prefix}__{ob.id}"
                segments[name] = int(bdf.loc[ob.id, "len"])
                seg_isolates[name].add(iso)
                named.append((name, ob.strand))

            for (n1, s1), (n2, s2) in zip(named, named[1:]):
                links.add((n1, s1, n2, s2))

    # Core anchors keep their true graph-wide occurrence depth (they are present
    # in every isolate regardless of which edges survive); accessory copies are
    # per-junction, so their depth is the isolate count within that junction.
    depths = {}
    for name, isos in seg_isolates.items():
        if "__" in name:
            depths[name] = len(isos)
        else:
            depths[name] = int(bdf.loc[name, "count"])
    return GFA(segments, links, depths), prefix_map
