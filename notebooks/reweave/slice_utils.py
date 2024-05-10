from classes import *


def slice_substitutions(i: Interval, S: Substitution):
    """Given a set of substitutions and an interval, returns the reduced set of
    substitutions relative to the interval."""
    new_S = []
    for s in filter(lambda s: i.position_is_in(s.pos), S):
        new_pos = s.pos - i.start
        new_S.append(Substitution(pos=new_pos, alt=s.alt))
    return new_S


def slice_deletions(i: Interval, D: Deletion):
    """Given a set of deletions and an interval, returns the reduced set of
    deletions relative to the interval."""
    new_D = []
    for d in filter(lambda d: i.overlap(d.pos, d.length), D):
        new_start = max(d.pos, i.start) - i.start
        new_end = min(d.pos + d.length, i.end) - i.start
        new_len = new_end - new_start
        new_D.append(Deletion(pos=new_start, length=new_len))
    return new_D


def slice_insertions(i: Interval, I: Insertion, block_L: int):
    """Given a set of insertions and an interval, returns the reduced set of
    insertions relative to the interval.
    Nb: insertions are left-inclusive, except for the last one, which is right-inclusive.
    """

    new_I = []
    for ins in filter(lambda ins: i.insertion_overlap(ins.pos, block_L), I):
        new_pos = ins.pos - i.start
        new_I.append(Insertion(pos=new_pos, ins=ins.ins))
    return new_I


def slice_edits(i: Interval, ed: Edit, block_L: int):
    """Given a set of edits and an interval, returns the reduced set of
    edits relative to the interval, together with the start/end coordinates
    of the interval on the node."""
    new_edits = Edit(
        ins=slice_insertions(i, ed.ins, block_L),
        dels=slice_deletions(i, ed.dels),
        subs=slice_substitutions(i, ed.subs),
    )
    return new_edits


def new_strandedness(old_strandedness, orientation, deep):
    """Given the strandedness of the previous node, the orientation of the match,
    and whether the block is the deep one, returns the strandedness of the new node.
    Inverted only if the block is shallow and the match orientation is inverted."""
    if deep or orientation:
        return old_strandedness
    else:
        return not old_strandedness


def new_position(
    old_position: tuple[int, int],
    node_coords: tuple[int, int],
    path_L: int,
    old_strandedness: bool,
):
    """
    This function takes as input:
    - the position (start, end) of the original node in the genome
    - the coordinates (ds, de) of the start/end of the interval projected on the node
    - the length of the genome
    - the strandedness of the original node
    And returns the new position of the node in the genome.
    Takes into account periodic boundary conditions for circular genomes.
    """
    old_s, old_e = old_position
    new_s_in_node, new_e_in_node = node_coords
    if old_strandedness:
        new_s = (old_s + new_s_in_node) % path_L
        new_e = (old_s + new_e_in_node) % path_L
    else:
        new_s = (old_e - new_e_in_node) % path_L
        new_e = (old_e - new_s_in_node) % path_L
    return new_s, new_e


def interval_node_coords(i: Interval, edits: Edit, block_L: int) -> tuple[int, int]:
    """Returns start/end coordinate of the interval on the node, including
    indels."""
    s, e = i.start, i.end
    for d in edits.dels:
        if d.pos <= i.start:
            l = min(d.length + d.pos, i.start) - d.pos
            s -= l
        if d.pos < i.end:
            l = min(d.length + d.pos, i.end) - d.pos
            e -= l
    for ins in edits.ins:
        if ins.pos < i.start:
            s += len(ins.ins)
        if ins.pos < i.end:
            e += len(ins.ins)
        # special case: insertion on the last position of the block and interval
        if (ins.pos == i.end) and (ins.pos == block_L):
            e += len(ins.ins)
    return s, e


def block_slice(b: Block, i: Interval, G: Pangraph) -> tuple[Block, dict[int, Node]]:
    """
    Given a block and an interval, slices the block along that interval.
    It returns:
    - a new block with the sliced consensus and alignment
    - a dictionary of {old node id -> new node}
    Nb: strandedness of new nodes is already ajdusted for post-merging, but the consensus
    sequence is always not changed. Reverse-complementation only occurs at merging.
    """

    # extract new consensus (always forward strand)
    new_consensus = b.consensus[i.start : i.end]
    block_L = b.consensus_len()

    node_updates = {}
    new_alignment = {}
    for old_node_id, edits in b.alignment.items():
        old_node = G.nodes[old_node_id]
        old_strandedness = old_node.strandedness
        new_strand = new_strandedness(old_strandedness, i.orientation, i.deep)
        path_L = G.paths[old_node.path_id].L
        node_coords = interval_node_coords(i, edits, block_L)
        new_pos = new_position(
            old_position=old_node.position,
            node_coords=node_coords,
            path_L=path_L,
            old_strandedness=old_strandedness,
        )

        # create new node
        new_node = Node(
            id=None,
            block_id=i.new_block_id,
            path_id=old_node.path_id,
            strandedness=new_strand,
            position=new_pos,
        )
        new_node.id = new_node.calculate_id()
        node_updates[old_node_id] = new_node

        # create new edits for alignment
        new_edits = slice_edits(i, edits, block_L)
        new_alignment[new_node.id] = new_edits

    # assemble new block
    new_block = Block(
        id=i.new_block_id, consensus=new_consensus, alignment=new_alignment
    )

    return new_block, node_updates
