from classes import *
from slice_utils import *
from intervals_utils import *
from collections import defaultdict
from itertools import groupby


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


def group_promises(H: list[tuple[Block, bool, bool]]) -> list[MergePromise]:
    """
    Given a list of blocks to merge, in the form of a list of tuples
    (block_slice, deep, orientation) it groups them by merge-id and
    returns a list of merge promises.
    """
    promises = []
    for new_bid, Bs in groupby(
        sorted(H, key=lambda x: x[0].id),
        key=lambda x: x[0].id,
    ):  # group by new block id
        # (in python I also need to sort first to ensure that all same ids come together)
        Bs = list(Bs)
        assert len(Bs) == 2, "Only two blocks can be merged"
        b1, deep1, strand1 = Bs[0]
        b2, deep2, strand2 = Bs[1]
        assert deep1 ^ deep2, "One block must be deep and the other shallow"
        assert strand1 == strand2, "Strandedness information is the same"
        b_deep = b1 if deep1 else b2
        b_shallow = b2 if deep1 else b1
        promises.append(
            MergePromise(b_deep=b_deep, b_shallow=b_shallow, orientation=strand1)
        )
    return promises
