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


def group_promises(H: list[ToMerge]) -> list[MergePromise]:
    """
    Given a list of blocks to merge, it groups them by block-id and
    returns a list of merge promises.
    Nb: the deep attribute is needed in case the two blocks have the same depth,
    then the tie is resolved by using the reference block as the deep one.
    This is done by the assign_deep_block function.
    """
    promises = []
    for new_bid, Bs in groupby(
        sorted(H, key=lambda x: x.block_id()),
        key=lambda x: x.block_id(),
    ):  # group by new block id
        # (in python I also need to sort first to ensure that all same ids come together)
        Bs = list(Bs)
        assert len(Bs) == 2, "Only two blocks can be merged"
        b1, b2 = Bs
        assert b1.deep ^ b2.deep, "One block must be deep and the other shallow"
        assert b1.orientation == b2.orientation, "orientation must be the same"
        b_deep = b1.block if b1.deep else b2.block
        b_shallow = b2.block if b1.deep else b1.block
        promises.append(
            MergePromise(b_deep=b_deep, b_shallow=b_shallow, orientation=b1.orientation)
        )
    return promises


def split_block(
    bid: int,
    M: list[Alignment],
    graph: Pangraph,
    thr_len: int,
) -> tuple[GraphUpdate, list[ToMerge]]:
    """
    This function takes as input:
    - a block id
    - a list of alignments
    - the pangenome graph
    - a threshold length for the blocks
    It takes case of splitting the block along the alignments, creating new blocks and nodes.
    blocks that need to be merged are included in a ToMerge object.
    The function also emits a GraphUpdate object, containing the instruction on how to
    update the graph after the block split.
    It returns:
    - a GraphUpdate object
    - a list of ToMerge objects
    """

    T = extract_hits(bid, M)

    # from hits to intervals
    L = graph.blocks[bid].consensus_len()
    I = extract_intervals(T, L, thr_len)

    # prepare graph update and list of blocks to be merged
    u = GraphUpdate(
        b_old_id=bid,
        b_new=[],
        n_new={nid: [] for nid in graph.blocks[bid].alignment.keys()},
    )
    H = []

    # cut block on intervals (generate new block/nodes)
    b = graph.blocks[bid]
    for i in I:
        b_slice, n_dict = block_slice(b, i, graph)

        # update graph update object
        u.b_new.append(b_slice)
        for old_nid, new_node in n_dict.items():
            u.n_new[old_nid].append(new_node)

        # add blocks that need to be merged to the list
        if i.aligned:
            H.append(
                ToMerge(
                    block=b_slice,
                    deep=i.deep,
                    orientation=i.orientation,
                )
            )

    # flip order in node lists for paths that are inverted
    for old_node_id, nodes in u.n_new.items():
        strand = graph.nodes[old_node_id].strandedness
        if not strand:
            u.n_new[old_node_id] = list(reversed(nodes))

    return u, H
