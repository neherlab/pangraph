from ..topology_utils import Node, Path, Edge


class JunctionNode(Node):
    """Node that also carries the unique node_id for unambiguous sequence lookup.

    Inherits equality/hashing from Node (compares block_id + strand only),
    so it works transparently with Path, Edge, and Junction.
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
        center: A Path of accessory blocks between the flanks.
        right: The core block on the right flank.
    """

    def __init__(self, left: Node, center: Path, right: Node) -> None:
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

    def flanks_bid(self, bid) -> bool:
        """Check whether the given block id is one of the flanking blocks."""
        return (self.left.id == bid) or (self.right.id == bid)

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

    def to_list(self):
        return [self.left.to_str_id(), self.center.to_list(), self.right.to_str_id()]

    @staticmethod
    def from_list(t) -> "Junction":
        return Junction(
            Node.from_str_id(t[0]), Path.from_list(t[1], circular=False), Node.from_str_id(t[2])
        )


def path_junction_split(path: Path, is_core) -> list[Junction]:
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
        path: A Path object (circular or linear).
        is_core: A callable taking a block id and returning True if the block is core.

    Returns:
        A list of Junction objects.
    """
    junctions = []

    current = []
    left_node = None
    for node in path.nodes:
        if is_core(node.id):
            J = Junction(left_node, Path(current), node)
            junctions.append(J)
            left_node = node
            current = []
        else:
            current.append(node)

    if path.circular:
        # complete periodic boundary: merge trailing non-core nodes into the first junction
        J = junctions[0]
        J.left = left_node
        J.center = Path(current + J.center.nodes)
        junctions[0] = J
    elif current or left_node is not None:
        # trailing accessory nodes after the last core block
        junctions.append(Junction(left_node, Path(current), None))

    return junctions
