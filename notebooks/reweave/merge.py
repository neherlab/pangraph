from classes import *
from utils import *


def solve_promise(p: MergePromise) -> Block:
    """
    Given a merge promise, it solves it by merging the two blocks.
    Returns the new block.
    """
    for n_id, edits in p.append_block.alignment.items():
        seq = apply_edits_to_ref(edits, p.appended_block.consensus)
        if p.orientation:
            # guarantees that the two sequences to be aligned have the same orientation
            seq = reverse_complement(seq)
        append_sequence(p.anchor_block, seq, n_id)
    return p.anchor_block
