from dataclasses import dataclass, field
from collections import defaultdict
from Bio.Seq import Seq
from Bio.Align import PairwiseAligner


def reverse_complement(seq: str) -> str:
    """Reverse complements a DNA sequence"""
    return str(Seq(seq).reverse_complement())


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
    circular: bool = None  # whether the genome is circular
    name: str = None  # name of the path


@dataclass
class Insertion:
    pos: int
    ins: str

    def reverse_complement(self, L) -> "Insertion":
        return Insertion(L - self.pos, reverse_complement(self.ins))

    def shift(self, shift: int) -> "Insertion":
        return Insertion(self.pos + shift, self.ins)


@dataclass
class Deletion:
    pos: int
    length: int

    def reverse_complement(self, L) -> "Deletion":
        return Deletion(L - self.pos - self.length, self.length)

    def shift(self, shift: int) -> "Deletion":
        return Deletion(self.pos + shift, self.length)


@dataclass
class Substitution:
    pos: int
    alt: str

    def reverse_complement(self, L) -> "Substitution":
        return Substitution(L - self.pos - 1, reverse_complement(self.alt))

    def shift(self, shift: int) -> "Substitution":
        return Substitution(self.pos + shift, self.alt)


@dataclass
class Edit:
    ins: list[Insertion] = field(default_factory=list)
    dels: list[Deletion] = field(default_factory=list)
    subs: list[Substitution] = field(default_factory=list)

    def reverse_complement(self, L) -> "Edit":
        return Edit(
            ins=[i.reverse_complement(L) for i in self.ins],
            dels=[d.reverse_complement(L) for d in self.dels],
            subs=[s.reverse_complement(L) for s in self.subs],
        )

    def shift(self, shift: int) -> "Edit":
        return Edit(
            ins=[i.shift(shift) for i in self.ins],
            dels=[d.shift(shift) for d in self.dels],
            subs=[s.shift(shift) for s in self.subs],
        )

    def concat(self, other: "Edit") -> "Edit":
        return Edit(
            ins=self.ins + other.ins,
            dels=self.dels + other.dels,
            subs=self.subs + other.subs,
        )


@dataclass
class Block:
    id: int  # block id
    consensus: str  # consensus sequence of the block
    alignment: dict[int, Edit]  # dict of node id to edit

    def depth(self) -> int:
        return len(self.alignment)

    def consensus_len(self) -> int:
        return len(self.consensus)

    def reverse_complement(self) -> "Block":
        """Returns a new block with the reverse complement of the consensus and alignment."""
        rev_cons = reverse_complement(self.consensus)
        L = self.consensus_len()
        rev_aln = {nid: e.reverse_complement(L) for nid, e in self.alignment.items()}
        return Block(self.id, rev_cons, rev_aln)


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
    anchor_block: str = None  # block used as anchor in a merger
    matches: int = None
    length: int = None
    quality: int = None
    cigar: str = None
    divergence: float = None
    align: float = None


@dataclass
class MergePromise:
    anchor_block: Block  # anchor block
    append_block: Block  # block that is appended to the anchor in the merge.
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

    def remove_path(self, pid: int):
        """Removes a path and all the corresponding nodes from the pangraph object."""
        nodes_to_remove = self.paths[pid].nodes
        for nid in nodes_to_remove:
            bid = self.nodes[nid].block_id
            # remove block alignment entry
            self.blocks[bid].alignment.pop(nid)
            # remove node from node dictionary
            self.nodes.pop(nid)
        # remove path from path dictionary
        self.paths.pop(pid)


@dataclass
class Interval:
    start: int
    end: int
    aligned: bool
    new_block_id: int
    is_anchor: bool = None
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
    is_anchor: bool
    orientation: bool

    def block_id(self) -> int:
        return self.block.id


def align_pairwise(ref: str, qry: str) -> tuple[str, str]:
    """
    Aligns two sequences using the global alignment algorithm.
    """
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 3
    aligner.mismatch_score = -1
    aligner.open_gap_score = -6
    aligner.extend_gap_score = 0
    aln = aligner.align(ref, qry)
    return aln[0]


def map_variations(consensus, qry):
    """
    Given a consensus and a query sequence, returns the edits required to transform the consensus into the query.
    Nb: re-implementation for testing purposes of:
    https://github.com/neherlab/pangraph/blob/d4272ece3ee326ca8228652667f9682f34d570a1/packages/pangraph/src/align/map_variations.rs#L8
    """
    aln = align_pairwise(consensus, qry)
    rc, qc = aln.coordinates
    edits = Edit()
    for i in range(len(rc) - 1):
        rb, qb = rc[i], qc[i]
        re, qe = rc[i + 1], qc[i + 1]
        if rb == re:
            # insertion
            edits.ins.append(Insertion(re, qry[qb:qe]))
        elif qb == qe:
            # deletion
            edits.dels.append(Deletion(rb, re - rb))
        else:
            # possible substitution
            r_seq = consensus[rb:re]
            q_seq = qry[qb:qe]
            if r_seq == q_seq:
                continue
            for j in range(len(r_seq)):
                if r_seq[j] != q_seq[j]:
                    edits.subs.append(Substitution(rb + j, q_seq[j]))
    return edits
