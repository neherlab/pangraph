"""Minimal, graph-agnostic GFA1 writer.

Given a set of segments (name -> length) and a set of links between them, write a
GFA1 file. Sequences are not stored: segment lengths are emitted as ``LN:i:``
tags and the sequence field is left as ``*``.
"""


def _orient(strand: bool) -> str:
    """Map a strand boolean to a GFA orientation symbol (True -> '+', False -> '-')."""
    return "+" if strand else "-"


def write_gfa(filepath, segments, links, depths=None) -> None:
    """Write a minimal GFA1 file.

    Args:
        filepath: Output path for the GFA file.
        segments: dict mapping segment name (str) -> length in bp (int).
        links: iterable of (from_name, from_strand, to_name, to_strand) tuples,
            where the strands are booleans (True -> '+', False -> '-'). Passing a
            set collapses duplicate links (e.g. the same edge seen on many
            isolates) automatically.
        depths: optional dict mapping segment name -> coverage depth. Emitted as a
            ``DP:f:`` tag (read by Bandage as node depth) for segments present in
            the dict.
    """
    depths = depths or {}
    with open(filepath, "w") as f:
        f.write("H\tVN:Z:1.0\n")
        for name, length in segments.items():
            line = f"S\t{name}\t*\tLN:i:{int(length)}"
            if name in depths:
                line += f"\tDP:f:{depths[name]}"
            f.write(line + "\n")
        for from_name, from_strand, to_name, to_strand in links:
            f.write(
                f"L\t{from_name}\t{_orient(from_strand)}"
                f"\t{to_name}\t{_orient(to_strand)}\t0M\n"
            )
