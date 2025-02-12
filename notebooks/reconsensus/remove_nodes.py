from utils import Pangraph

# remove empty nodes:
# - for blocks that originate from a new merger or split
# - goes through each block and removes any node with an empty sequence
# - also remove the node from the graph node list and paths


def find_empty_nodes(graph: Pangraph, block_ids: list[int]) -> list[int]:
    """
    Finds nodes with empty sequences (full deletion) in the graph.
    It only checks specific blocks (the ones that were updated by a merger/split).
    """
    node_ids_to_delete = []
    for block_id in block_ids:
        block = graph.blocks[block_id]
        L = len(block.consensus)
        for node_id, edits in block.alignment.items():
            if len(edits.ins) != 0 or len(edits.subs) != 0 or len(edits.dels) == 0:
                continue
            if sum([d.length for d in edits.dels]) == L:
                node_ids_to_delete.append(node_id)
                # debug assert: node should have size zero
                node = graph.nodes[node_id]
                assert node.position[0] == node.position[1], "node has non-zero size"

    return node_ids_to_delete


def remove_nodes_from_graph(graph: Pangraph, node_ids: list[int]):
    """
    Removes each node from the graph inplace. The node gets removed from:
    - the node dictionary
    - the alignment object in block
    - the path
    """
    for node_id in node_ids:
        node = graph.nodes[node_id]
        path_id = node.path_id
        block_id = node.block_id

        # remove from node dictionary
        del graph.nodes[node_id]

        # remove from path
        path_nodes = graph.paths[path_id].nodes
        node_idx = path_nodes.index(node_id)
        path_nodes.pop(node_idx)

        # remove from block alignment
        del graph.blocks[block_id].alignment[node_id]


def remove_emtpy_nodes(graph: Pangraph, block_ids: list[int]):
    """
    If present, removes empty nodes from a graph inplace. Only specified
    blocks are checked for empty nodes.
    """
    node_ids = find_empty_nodes(graph, block_ids)
    remove_nodes_from_graph(graph, node_ids)
