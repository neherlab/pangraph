from dataclasses import dataclass, field
from collections import defaultdict


@dataclass
class Node:
    id: int  # node id
    block_id: int  # block id
    path_id: int  # path id
    position: tuple[int]  # start/end position on the genome
    strandedness: bool  # strandedness of the node

    def calculate_id(self) -> int:
        """the id of the node is given by the hash of the tuple (block_id, path_id, start, end)"""
        nid = hash((self.block_id, self.path_id, self.position[0], self.position[1]))
        return nid


@dataclass
class Path:
    id: int  # path id
    nodes: list[int]  # list of node ids
    L: int  # total length of the genome


@dataclass
class Insertion:
    pos: int
    ins: str


@dataclass
class Deletion:
    pos: int
    length: int


@dataclass
class Substitution:
    pos: int
    alt: str


@dataclass
class Edit:
    ins: list[Insertion] = field(default_factory=list)
    dels: list[Deletion] = field(default_factory=list)
    subs: list[Substitution] = field(default_factory=list)


def apply_edits_to_ref(edits: Edit, ref: str) -> str:
    """
    Apply the edits to the reference to obtain the query sequence
    """
    qry = list(ref)
    for S in edits.subs:
        qry[S.pos] = S.alt
    for D in edits.dels:
        for l in range(D.length):
            qry[D.pos + l] = ""
    for I in edits.ins:
        if I.pos > 0:
            qry[I.pos - 1] += I.ins
        elif I.pos == 0:
            qry[0] = I.ins + qry[0]
    return "".join(qry)


def reverse_complement(seq: str) -> str:
    """
    Reverse complement a DNA sequence
    """
    complement = {"A": "T", "C": "G", "G": "C", "T": "A"}
    return "".join(complement[base] for base in reversed(seq))


@dataclass
class Block:
    id: int  # block id
    consensus: str  # consensus sequence of the block
    alignment: dict[int, Edit]  # dict of node id to edit

    def depth(self) -> int:
        return len(self.alignment)

    def consensus_len(self) -> int:
        return len(self.consensus)


@dataclass
class Hit:
    name: int  # name/id of the block
    length: int  # length of the block (not the hit!)
    start: int  # his start position
    stop: int  # hit stop position (exclusive)


@dataclass
class Alignment:
    qry: Hit
    reff: Hit
    orientation: bool  # fwd/rev match orientation
    new_block_id: int = None  # id of the newly-created block. Also alignment id.
    deep_block: str = None  # decide whether qry or ref is the deepest block.
    matches: int = None
    length: int = None
    quality: int = None
    cigar: str = None
    divergence: float = None
    align: float = None


@dataclass
class MergePromise:
    b_deep: Block  # deep block
    b_shallow: Block  # shallow block (can be reverse-complemented)
    orientation: bool  # fwd/rev match orientation


@dataclass
class GraphUpdate:
    b_old_id: int
    b_new: list[Block]
    n_new: dict[int, list[Node]]
    # nb: node list is already in the order of the new path


@dataclass
class Pangraph:
    paths: dict[int, Path]  # dict of path id to path
    blocks: dict[int, Block]  # dict of block id to block
    nodes: dict[int, Node]  # dict of node id to node

    def update(self, u: GraphUpdate):
        """Applies a graph update to the pangraph object.
        Nb: the updated already contains the correct order of nodes."""

        # consistency check: node ids
        old_nodes_set_from_graph = set(self.blocks[u.b_old_id].alignment.keys())
        old_nodes_set_from_update = set(u.n_new.keys())
        assert (
            old_nodes_set_from_graph == old_nodes_set_from_update
        ), f"old nodes mismatch: {old_nodes_set_from_graph} != {old_nodes_set_from_update}"

        # remove old block and add new ones
        self.blocks.pop(u.b_old_id)
        self.blocks.update({b.id: b for b in u.b_new})

        for old_nid, new_nodes in u.n_new.items():
            # find path
            path_id = self.nodes[old_nid].path_id

            # remove old node from path
            old_idx = self.paths[path_id].nodes.index(old_nid)
            self.paths[path_id].nodes.pop(old_idx)
            # add new nodes to path
            new_ids = [n.id for n in new_nodes]
            self.paths[path_id].nodes[old_idx:old_idx] = new_ids

            # remove old node from graph node dictionary
            self.nodes.pop(old_nid)

            # add new nodes to graph node dictionary
            self.nodes.update({n.id: n for n in new_nodes})

    # def get_node_seq(self, nid):
    #     """Get the sequence of a node.
    #     The strandedness is the one of the block, not of the path.
    #     It might be necessary to reverse-complement it to get the
    #     path sequence."""
    #     bid = self.nodes[nid].block_id
    #     b = self.blocks[bid]
    #     edit = b.alignment[nid]
    #     return apply_edits_to_ref(edit, b.consensus)

    # def get_path_seq(self, pid):
    #     """Reconstruct the full sequence of a path."""
    #     path = self.paths[pid]
    #     seq = ""
    #     for nid in path.nodes:
    #         n = self.nodes[nid]
    #         strand = n.strandedness
    #         node_seq = self.get_node_seq(nid)
    #         if not strand:
    #             node_seq = reverse_complement(node_seq)
    #         seq += node_seq
    #     return seq


@dataclass
class Interval:
    start: int
    end: int
    aligned: bool
    new_block_id: int
    deep: bool = None
    orientation: bool = None

    def __len__(self) -> int:
        assert (
            self.end >= self.start
        ), f"end position {self.end} < start position {self.start}"
        return self.end - self.start

    def position_is_in(self, pos: int) -> bool:
        """Checks whether the position is on the interval."""
        return self.start <= pos < self.end

    def overlap(self, s: int, l: int) -> bool:
        """Given the start and length of an interval, returns whether there is an overlap."""
        if self.position_is_in(s):
            return True
        if s < self.start:
            if s + l > self.start:
                return True
            else:
                return False
        return False

    def insertion_overlap(self, ins_pos: int, block_L: int) -> bool:
        """Decides whether an insertion overlaps with the interval.
        insertions are left-inclusive with one exception:
        the one at the righ-edge of the block is included if self.end == block_L.
        """
        if self.position_is_in(ins_pos):
            return True
        elif (ins_pos == block_L) and (self.end == block_L):
            return True
        else:
            return False


@dataclass
class ToMerge:
    block: Block
    deep: bool
    orientation: bool

    def block_id(self) -> int:
        return self.block.id
