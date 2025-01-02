from dataclasses import dataclass


@dataclass
class SimpleNode:
    bid: int  # block id
    strand: bool  # strand

    def __invert__(self) -> "SimpleNode":
        return SimpleNode(self.bid, not self.strand)

    def __eq__(self, other: object) -> bool:
        return self.bid == other.bid and self.strand == other.strand

    def __hash__(self) -> int:
        return hash((self.bid, self.strand))

    def __repr__(self) -> str:
        s = "+" if self.strand else "-"
        return f"({self.bid}{s})"


# already implemented in
# https://github.com/neherlab/pangraph/blob/98886771cb20cd4bfe7ce33c52dafc2fc33f6faa/packages/pangraph/src/circularize/circularize_utils.rs#L41
@dataclass
class Edge:
    n1: SimpleNode
    n2: SimpleNode

    def __invert__(self) -> "Edge":
        return Edge(~self.n2, ~self.n1)

    def oriented_equal(self, other: "Edge") -> bool:
        return self.n1 == other.n1 and self.n2 == other.n2

    def __eq__(self, other: object) -> bool:
        return self.oriented_equal(other) or self.oriented_equal(~other)

    def side_hash(self) -> int:
        return hash((self.n1, self.n2))

    def __hash__(self) -> int:
        return self.side_hash() ^ (~self).side_hash()

    def __repr__(self) -> str:
        return f"[{self.n1}|{self.n2}]"


@dataclass
class SimplePath:
    nodes: list[SimpleNode]
    circular: bool

    def to_edges(self) -> list[Edge]:
        edges = [
            Edge(self.nodes[i], self.nodes[i + 1]) for i in range(len(self.nodes) - 1)
        ]
        if self.circular:
            edges.append(Edge(self.nodes[-1], self.nodes[0]))
        return edges
