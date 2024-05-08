from classes import *
from collections import defaultdict


def assign_new_block_ids(mergers: list[Alignment]):
    """assigns new block id to each merger inplace"""
    for i, a in enumerate(mergers):
        # this could also be done with a sequential counter if we do not parallelize the tree traversal
        q = a.qry
        r = a.reff
        vec = (q.name, q.start, q.stop, r.name, r.start, r.stop, i)
        a.new_block_id = hash(vec)


def assign_deep_block(mergers: list[Alignment], graph: Pangraph):
    """
    Decide whether ref or qry is the deep block.
    We will append sequences to this block, and it will not be flipped if
    the match is on the reverse strand.
    """
    for m in mergers:
        ref_block = graph.blocks[m.reff.name]
        qry_block = graph.blocks[m.qry.name]
        ref_depth = ref_block.depth()
        qry_depth = qry_block.depth()
        if ref_depth >= qry_depth:
            m.deep_block = "reff"
        else:
            m.deep_block = "qry"


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


def extract_hits(bid: int, M: list[Alignment]) -> list[dict]:
    """Given a block id and list of alignments, returns a list of hits
    on that block, with information on strandedness and whether the block is
    the deepest one of the pair."""
    hits = []
    for m in M:
        if m.reff.name == bid:
            deep = m.deep_block == "reff"
            hits.append(
                {
                    "new_block_id": m.new_block_id,
                    "deep": deep,
                    "orientation": m.orientation,
                    "hit": m.reff,
                }
            )
        if m.qry.name == bid:
            deep = m.deep_block == "qry"
            hits.append(
                {
                    "new_block_id": m.new_block_id,
                    "deep": deep,
                    "orientation": m.orientation,
                    "hit": m.qry,
                }
            )
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
        deep=hit_dict["deep"],
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
