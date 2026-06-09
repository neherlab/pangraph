"""Minimal, graph-agnostic GFA1 representation.

An in-memory GFA1 graph (segments + links, with optional per-segment depths)
that knows how to serialize itself. Sequences are not stored: segment lengths are
emitted as ``LN:i:`` tags and the sequence field is left as ``*``. The real bp
lengths can optionally be rescaled at write time (e.g. ``log`` or ``1/l``) for
visualization — see :meth:`GFA.write`.
"""


def _orient(strand: bool) -> str:
    """Map a strand boolean to a GFA orientation symbol (True -> '+', False -> '-')."""
    return "+" if strand else "-"


class GFA:
    """An in-memory GFA1 graph: segments, links, and optional per-segment depths.

    Build it directly or via a builder (e.g.
    :func:`pypangraph.export.junction_context_gfa`), then call :meth:`write` to
    serialize it to a GFA1 file.
    """

    def __init__(self, segments, links, depths=None) -> None:
        """Store the graph.

        Args:
            segments: dict mapping segment name (str) -> length in bp (int).
            links: iterable of (from_name, from_strand, to_name, to_strand) tuples,
                where the strands are booleans (True -> '+', False -> '-'). A set
                collapses duplicate links automatically.
            depths: optional dict mapping segment name -> coverage depth, emitted
                as a ``DP:f:`` tag (read by Bandage as node depth).
        """
        self.segments = segments
        self.links = links
        self.depths = depths or {}

    def write(self, filepath, length_transform=None) -> None:
        """Write this graph to a minimal GFA1 file.

        Emits the ``H`` header, one ``S`` line per segment (``LN:i:`` length and,
        when present in :attr:`depths`, a ``DP:f:`` depth tag, sequence ``*``),
        and one ``L`` line per link with a ``0M`` overlap.

        Args:
            filepath: Output path for the GFA1 file.
            length_transform: optional ``Callable`` applied to each segment's real
                bp length to rescale the emitted ``LN:i:`` value (useful for
                visualization, since block lengths span orders of magnitude). The
                result is rounded to the nearest integer and clamped to a minimum
                of ``1``; ``None`` (default) emits the true length. The caller owns
                the scaling: raw ``math.log(l)`` is tiny (~5-14) and raw ``1/l``
                collapses everything to the ``1`` floor, so scale as needed, e.g.::

                    gfa.write("out.gfa", length_transform=lambda l: 100 * math.log(l))
                    gfa.write("out.gfa", length_transform=lambda l: 1e7 / l)
        """
        with open(filepath, "w") as f:
            f.write("H\tVN:Z:1.0\n")
            for name, length in self.segments.items():
                ln = (
                    int(length)
                    if length_transform is None
                    else max(1, round(length_transform(length)))
                )
                line = f"S\t{name}\t*\tLN:i:{ln}"
                if name in self.depths:
                    line += f"\tDP:f:{self.depths[name]}"
                f.write(line + "\n")
            for from_name, from_strand, to_name, to_strand in self.links:
                f.write(
                    f"L\t{from_name}\t{_orient(from_strand)}"
                    f"\t{to_name}\t{_orient(to_strand)}\t0M\n"
                )
