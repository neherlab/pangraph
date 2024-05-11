from classes import *


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


def __unaligned_interval(start: int, end: int, block_id: int):
    return Interval(
        start=start,
        end=end,
        aligned=False,
        new_block_id=hash((block_id, start, end)),
    )


def __aligned_interval(hit_dict: dict):
    return Interval(
        start=hit_dict["hit"].start,
        end=hit_dict["hit"].stop,
        aligned=True,
        new_block_id=hit_dict["new_block_id"],
        is_anchor=hit_dict["is_anchor"],
        orientation=hit_dict["orientation"],
    )


def create_intervals(
    hits: list[dict],
    block_L: int,
) -> list[Interval]:
    """Create intervals from the list of hits. This split the block in alternated
    aligned and unaligned intervals.
    Assigns new block ids to non-merged intervals."""

    # sort hits by start position
    hits.sort(key=lambda x: x["hit"].start)

    I = []  # list of intervals
    cursor = 0
    for h_dict in hits:
        h = h_dict["hit"]
        if h.start > cursor:
            # unaligned interval
            I.append(__unaligned_interval(start=cursor, end=h.start, block_id=h.name))
            cursor = h.start
        # aligned interval
        I.append(__aligned_interval(h_dict))
        cursor = h.stop
    if cursor < block_L:
        # final unaligned interval
        I.append(__unaligned_interval(start=cursor, end=block_L, block_id=h.name))

    return I


def refine_intervals(I: list[Interval], thr_len: int) -> list[Interval]:
    """Intervals shorter than the threshold length are joined to the longest adjacent aligned interval."""
    mergers = []
    for n, i in enumerate(I):
        # unaligned intervals shorter than the threshold are joined to the longest adjacent aligned interval
        if len(i) < thr_len:

            # find the longest adjacent aligned interval
            L_left = len(I[n - 1]) if n > 0 else 0
            L_right = len(I[n + 1]) if n + 1 < len(I) else 0

            # sanity-check: aligned interval should always be longer than the threshold.
            assert not i.aligned, f"aligned interval {n} shorter than the threshold len"
            # sanity check: left and right intervals should be aligned and longer than threshold
            if n > 0:
                assert I[n - 1].aligned, f"no adjacent aligned interval on the left"
                assert L_left >= thr_len, f"left interval shorter than threshold len"
            if n + 1 < len(I):
                assert I[n + 1].aligned, f"no adjacent aligned interval on the right"
                assert L_right >= thr_len, f"right interval shorter than threshold len"

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


def extract_intervals(
    hits: list[tuple[tuple[int, str], Hit]], block_L: int, thr_len: int
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
