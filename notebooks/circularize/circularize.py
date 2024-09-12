from collections import defaultdict, Counter
from dataclasses import dataclass
from utils import Node, Path, Pangraph, Block


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


def block_depths(graph: Pangraph) -> dict[int, int]:
    """
    Returns a dictionary block id -> block depth.
    """
    return {bid: len(b.alignment) for bid, b in graph.blocks.items()}


def graph_to_simplepaths(graph: Pangraph) -> dict[int, SimplePath]:
    """
    Returns a dictionary path_id -> SimplePath.
    """
    simple_paths = {}
    for path_id, path in graph.paths.items():
        simple_nodes = []
        for node_id in path.nodes:
            node = graph.nodes[node_id]
            simple_node = SimpleNode(node.block_id, node.strandedness)
            simple_nodes.append(simple_node)
        simple_path = SimplePath(simple_nodes, path.circular)
        simple_paths[path_id] = simple_path
    return simple_paths


def edges_count(graph: Pangraph) -> dict[Edge, int]:
    """
    Returns a dictionary of edge -> edge count.
    """
    paths = graph_to_simplepaths(graph)
    edges_ct = Counter()
    for path in paths.values():
        edges_ct.update(path.to_edges())
    return dict(edges_ct)


def transitive_edges(graph: Pangraph) -> list[Edge]:
    """
    Returns a list of transitive edges in the graph.
    """
    n_blocks = block_depths(graph)
    n_edges = edges_count(graph)
    transitive_edges = []
    for edge, ct in n_edges.items():
        bid1 = edge.n1.bid
        bid2 = edge.n2.bid
        if n_blocks[bid1] == n_blocks[bid2] == ct:
            transitive_edges.append(edge)
    return transitive_edges
