from collections import Counter
from utils import Pangraph
from merge_blocks import merge_blocks
from circularize_utils import Edge, SimpleNode, SimplePath


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
    Returns a list of transitive edges between two different blocks in the graph
    (no self-loops).
    """
    n_blocks = block_depths(graph)
    n_edges = edges_count(graph)
    transitive_edges = []
    for edge, ct in n_edges.items():
        bid1 = edge.n1.bid
        bid2 = edge.n2.bid
        if n_blocks[bid1] == n_blocks[bid2] == ct and bid1 != bid2:
            transitive_edges.append(edge)
    return transitive_edges


def remove_transitive_edges(graph: Pangraph) -> Pangraph:
    """
    Removes transitive edges from the graph inplace.
    """
    # TODO: can be improved: if transitive edges are between separate nodes,
    # no need to recalculate the list of transitive edges every time
    while len(t_edges := transitive_edges(graph)) > 0:
        merge_blocks(graph, t_edges[0])
