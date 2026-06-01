"""Linear schematic plot of a junction across isolates."""

import random
from collections import defaultdict


def _random_color() -> tuple[float, float, float]:
    """Return a mid-saturation RGB triple, each channel uniform in [0.3, 0.95]."""
    return (
        random.uniform(0.3, 0.95),
        random.uniform(0.3, 0.95),
        random.uniform(0.3, 0.95),
    )


def linear_junction_plot(
    ax,
    bj,
    edge: str,
    *,
    isolates: list[str] | None = None,
    color_map: dict | None = None,
    left_flank_color="C0",
    right_flank_color="C1",
    highlight_inverted: bool = False,
) -> dict:
    """Draw a per-isolate linear schematic of one core-edge junction.

    Each isolate becomes a row of horizontal bars (one per oriented block),
    with bar width set to the block's consensus length. Junctions are
    co-oriented via ``Junction.to_canonical()`` so the left flank lines up
    across rows.

    Args:
        ax: A Matplotlib axes object; mutated in place.
        bj: A ``BackboneJunctions`` instance.
        edge: Canonical edge string ID (e.g. ``"100_f__200_f"``).
        isolates: Optional explicit row order. Defaults to
            ``sorted(bj[edge].keys())``. Unknown names raise ``KeyError``.
        color_map: Optional seed mapping ``block_id -> color`` for accessory
            blocks. Missing entries are filled in with random mid-saturated
            colors. The input is not mutated. Flank colors always override
            this mapping.
        left_flank_color: Face color of the left-flank block.
        right_flank_color: Face color of the right-flank block.
        highlight_inverted: When ``True``, blocks with negative strand get
            a red border; otherwise all borders are black.

    Returns:
        A plain ``dict`` mapping accessory block id to color, including any
        newly-generated random colors. Pass it as ``color_map`` to a
        subsequent call to keep block colors consistent across panels.
    """
    per_iso = bj[edge]

    if isolates is None:
        isolates = sorted(per_iso.keys())

    seed = color_map or {}
    colors = defaultdict(_random_color, seed)

    block_len = bj._bdf["len"]

    for row, iso in enumerate(isolates):
        J = per_iso[iso].to_canonical()
        obs = J.oriented_blocks()
        n = len(obs)
        x = 0
        for i, ob in enumerate(obs):
            length = block_len.loc[ob.id]
            if i == 0:
                face = left_flank_color
            elif i == n - 1:
                face = right_flank_color
            else:
                face = colors[ob.id]
            edge_c = "red" if (highlight_inverted and not ob.strand) else "black"
            ax.barh(
                row,
                length,
                left=x,
                height=0.8,
                color=face,
                edgecolor=edge_c,
                linewidth=0.4,
            )
            x += length

    ax.set_yticks(range(len(isolates)))
    ax.set_yticklabels(isolates, fontsize=8)
    ax.invert_yaxis()
    ax.set_xlabel("position along junction (bp)")

    return dict(colors)
