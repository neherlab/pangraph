from __future__ import annotations

from ..topology_utils import OrientedBlock, Walk, Edge


class JunctionNode(OrientedBlock):
    """OrientedBlock that also carries the unique node_id for unambiguous sequence lookup.

    Inherits equality/hashing from OrientedBlock (compares block_id + strand only),
    so it works transparently with Walk, Edge, and Junction.
    """

    def __init__(self, bid, strand: bool, node_id: int) -> None:
        super().__init__(bid, strand)
        self.node_id = node_id

    def invert(self) -> "JunctionNode":
        return JunctionNode(self.id, not self.strand, self.node_id)

    def __repr__(self) -> str:
        s = "+" if self.strand else "-"
        return f"[{self.id}|{s}|n{self.node_id}]"


class Junction:
    """A junction is a segment of accessory blocks flanked by two core blocks,
    with reverse-complement symmetry.

    Attributes:
        left: The core block on the left flank.
        center: A Walk of accessory blocks between the flanks.
        right: The core block on the right flank.
    """

    def __init__(self, left: OrientedBlock, center: Walk, right: OrientedBlock) -> None:
        self.left = left
        self.center = center
        self.right = right

    def invert(self) -> "Junction":
        return Junction(self.right.invert(), self.center.invert(), self.left.invert())

    def flanking_edge(self) -> Edge | None:
        """Returns the Edge connecting the two flanking core blocks.

        Returns None for terminal junctions on linear paths where one flank is missing.
        """
        if self.left is None or self.right is None:
            return None
        return Edge(self.left, self.right)

    def is_canonical(self, edge: Edge) -> bool:
        """Whether this junction is oriented along the edge's canonical direction.

        A junction for a given edge can appear in either orientation across genomes;
        it is canonical when its left flank matches the edge's left node. Returns False
        for terminal junctions with no left flank.
        """
        if self.left is None:
            return False
        return self.left == edge.left

    def __side_eq__(self, o: object) -> bool:
        return self.left == o.left and self.center == o.center and self.right == o.right

    def __eq__(self, o: object) -> bool:
        return self.__side_eq__(o) or self.__side_eq__(o.invert())

    def __side_hash__(self) -> int:
        return hash((self.left, self.center, self.right))

    def __hash__(self) -> int:
        return self.__side_hash__() ^ self.invert().__side_hash__()

    def __repr__(self) -> str:
        return f"{self.left} <-- {self.center} --> {self.right}"


def path_junction_split(path: Walk, is_core) -> list[Junction]:
    """Split a path into junctions at core block boundaries.

    Given a path and a boolean predicate on node ids, splits the path into a list
    of Junctions. Each junction has flanking core blocks (for which the predicate
    returns True) and a center path of non-core blocks.

    For circular paths, the last junction wraps around the origin, connecting the
    last core block back to the first.

    For linear paths, there may be terminal junctions at the start and/or end with
    only one flanking core block (left=None or right=None). These terminal junctions
    have no flanking edge.

    Args:
        path: A Walk object (circular or linear).
        is_core: A callable taking a block id and returning True if the block is core.

    Returns:
        A list of Junction objects.

    Raises:
        ValueError: if the path has fewer than two core blocks, in which case no
            junction can be defined (and a circular path would have no flanks at all).
    """
    n_core = sum(1 for node in path.nodes if is_core(node.id))
    if n_core < 2:
        raise ValueError(
            f"path has {n_core} core block(s); at least 2 are required to define a junction"
        )

    # collect here the junctions
    junctions = []

    current = []
    left_node = None
    for node in path.nodes:
        if is_core(node.id):
            J = Junction(left_node, Walk(current), node)
            junctions.append(J)
            left_node = node
            current = []
        else:
            current.append(node)

    if path.circular:
        # complete periodic boundary: merge trailing non-core nodes into the first junction
        J = junctions[0]
        J.left = left_node
        J.center = Walk(current + J.center.nodes)
        junctions[0] = J
    elif current or left_node is not None:
        # trailing accessory nodes after the last core block
        junctions.append(Junction(left_node, Walk(current), None))

    return junctions
