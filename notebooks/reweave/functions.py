from classes import *
from collections import defaultdict


def target_blocks(mergers: list[Alignment]) -> dict[int, list[Alignment]]:
    """
    given a list of mergers, returns a dictionary of target blocks to
    the alignments that target them.
    """
    target_blocks = defaultdict(list)
    for merger in mergers:
        qry = merger.qry.name
        reff = merger.reff.name
        target_blocks[qry].append(merger)
        target_blocks[reff].append(merger)
    return target_blocks


def extract_hits(bid: int, M: list[Alignment]) -> list[tuple[int, Hit]]:
    """Collect all of the hits on the block, and returns a list of tuples (id, hit) for each hit.
    Hits are returned sorted by start position.
    Nb: can be more than a single hit per alignment, in case of a self-map."""
    hits = []
    for m in M:
        if m.qry.name == bid:
            hits.append((m.id, m.qry, "qry"))
        if m.reff.name == bid:
            hits.append((m.id, m.reff, "reff"))
    return hits


def intervals_sanity_checks(I: list[Interval], L: int) -> None:
    """
    Sanity checks to verify that:
    - intervals cover the whole block
    - intervals are contiguous
    - two consecutive intervals are not both unaligned.
    """
    assert I[0].start == 0, f"first interval does not start at 0"
    assert I[-1].end == L, f"last interval does not end at the block length"
    for n in range(1, len(I)):
        assert I[n - 1].end == I[n].start, f"interval {n-1} and {n} are not contiguous"
        assert (
            I[n - 1].aligned or I[n].aligned
        ), f"two consecutive unaligned intervals: {n-1} and {n}"


def create_intervals(hits: list[tuple[int, Hit]], block_L: int) -> list[Interval]:
    """Create intervals from the list of hits. This split the block in alternated
    aligned and unaligned intervals."""

    # sort hits by start position
    hits.sort(key=lambda x: x[1].start)

    I = []  # list of intervals
    cursor = 0
    for match_id, h in hits:
        if h.start > cursor:
            # unaligned interval
            I.append(Interval(cursor, h.start, False))
            cursor = h.start
        # aligned interval
        I.append(Interval(h.start, h.stop, True, match_id))
        cursor = h.stop
    if cursor < block_L:
        # final unaligned interval
        I.append(Interval(cursor, block_L, False))
    return I


def extract_intervals(
    hits: list[tuple[int, Hit]], block_L: int, thr_len: int
) -> list[Interval]:
    """Split the block into intervals, according to the hits.
    block_L is the total block consensus length."""
    # create intervals
    intervals = create_intervals(hits, block_L)

    # refine intervals inplace: merges short unaligned intervals with the longest adjacent aligned interval
    refine_intervals(intervals, thr_len)

    # assertions for sanity checks
    intervals_sanity_checks(intervals, block_L)

    return intervals


def refine_intervals(I: list[Interval], thr_len: int) -> list[Interval]:
    """Intervals shorter than the threshold length are joined to the longest adjacent aligned interval."""
    mergers = []
    for n, i in enumerate(I):
        # unaligned intervals shorter than the threshold are joined to the longest adjacent aligned interval
        if len(i) < thr_len:

            # sanity-check: aligned interval should always be longer than the threshold.
            assert (
                not i.aligned
            ), f"aligned interval {n} is shorter than the threshold length"

            # find the longest adjacent aligned interval
            L_left = len(I[n - 1]) if n > 0 else 0
            L_right = len(I[n + 1]) if n + 1 < len(I) else 0
            if n > 0:
                assert I[n - 1].aligned, f"no adjacent aligned interval on the left"
            if n + 1 < len(I):
                assert I[n + 1].aligned, f"no adjacent aligned interval on the right"
            assert (L_left > thr_len) or (
                L_right > thr_len
            ), f"no adjacent aligned interval longer than the threshold length"

            # merge with the longest adjacent aligned interval
            # if same length: merge with the left interval
            mergers.append((n, n - 1) if L_left >= L_right else (n, n + 1))

    # merge intervals. To avoid index issues, we start from the end of the list.
    # and we first modify the adjacent interval, and then remove the current one.
    for n_from, n_to in mergers[::-1]:
        if n_from < n_to:
            I[n_to].start = I[n_from].start
        else:
            I[n_to].end = I[n_from].end
        del I[n_from]
