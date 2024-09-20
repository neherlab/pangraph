from .utils import Pangraph, Block, reverse_complement
from .circularize_utils import Edge, SimpleNode
from .utils import Node
from notebooks.dump import dump


def merge_blocks(graph: Pangraph, edge: Edge):
    """
    Merges the two blocks from the same edge into a single block,
    modifying the graph in place.
    """

    # 0) orient edge: deterministically assign a left and right block
    edge = orient_merging_edge(graph, edge)

    # 1) find node pairs: recovers pairs of nodes to merge,
    # and creates new nodes to be added to the graph
    node_map, new_nodes = find_node_pairings(graph, edge)

    # 2) merge alignment: concatenates the alignment of the two blocks
    # and returns a new block
    new_block = merge_alignment(graph, edge, node_map, new_nodes)

    # 3) update paths and nodes
    graph_merging_update(graph, new_block, new_nodes, edge)


def orient_merging_edge(graph: Pangraph, edge: Edge) -> Edge:
    """
    Returns the oriented edge for merging, with the first block being the anchor block.
    This is:
    - the block with the longest consensus sequence
    - if the consensus lengths are equal, the block with the smallest id
    """
    b1 = graph.blocks[edge.n1.bid]
    b2 = graph.blocks[edge.n2.bid]
    L1, L2 = b1.consensus_len(), b2.consensus_len()
    if L1 > L2:
        return edge
    elif L1 < L2:
        return ~edge
    elif b1.id < b2.id:
        return edge
    else:
        return ~edge


def find_node_pairings(graph: Pangraph, edge: Edge) -> dict[int, int]:
    """
    Given a specific edge between blocks, returns pairings for all nodes in the blocks,
    and a dictionary of new nodes to be added to the graph.
    """
    node_pairings = {}
    new_nodes = {}
    for path_id, path in graph.paths.items():
        N = len(path.nodes)  # number of nodes
        I = N if path.circular else N - 1  # number of edges
        for i in range(I):
            # get node ids, nodes, block ids and strandedness of the two nodes
            nid1, nid2 = path.nodes[i], path.nodes[(i + 1) % N]
            n1, n2 = graph.nodes[nid1], graph.nodes[nid2]
            bid1, bid2 = n1.block_id, n2.block_id
            strand1, strand2 = n1.strandedness, n2.strandedness

            # check if the two nodes correspond to the desired edge
            sn1, sn2 = SimpleNode(bid1, strand1), SimpleNode(bid2, strand2)
            e = Edge(sn1, sn2)
            if edge == e:
                # create node pairings
                node_pairings[nid1] = nid2
                node_pairings[nid2] = nid1

                # create new node, in the orientation of the left edge
                s1, e1 = n1.position
                s2, e2 = n2.position
                if edge.n1 == sn1:
                    # edge orientation is the same as the path
                    new_s = s1
                    new_e = e2
                    assert (e1 % path.L) == (s2 % path.L), "nodes must be adjacent"
                    new_strand = strand1
                elif edge.n1 == ~sn2:
                    # edge orientation is the opposite of the path
                    new_s = s2
                    new_e = e1
                    assert (e2 % path.L) == (s1 % path.L), "nodes must be adjacent"
                    new_strand = strand2
                else:
                    raise ValueError("unexpected edge orientation")
                # create new node and assign new id based on hash
                new_node = Node(None, edge.n1.bid, path_id, (new_s, new_e), new_strand)
                new_node.id = new_node.calculate_id()
                new_nodes[nid1] = new_node
                new_nodes[nid2] = new_node

    return node_pairings, new_nodes


def concatenate_alignments(bl1, bl2, node_map, new_nodes_id, new_block_id):
    """Concatenates two blocks, with new node ids."""
    seq1, seq2 = bl1.consensus, bl2.consensus
    aln1, aln2 = bl1.alignment, bl2.alignment
    L1 = bl1.consensus_len()

    assert bl1.depth() == bl2.depth(), "blocks must have the same depth"

    seq = seq1 + seq2
    aln = {}
    for nid1, e1 in aln1.items():
        nid2 = node_map[nid1]
        e2 = aln2[nid2]
        new_id = new_nodes_id[nid1]
        aln[new_id] = e1.concat(e2.shift(L1))
    return Block(new_block_id, seq, aln)


def merge_alignment(
    graph: Pangraph, edge: Edge, node_map: dict[int, int], new_nodes: dict[int, Node]
):
    """
    Merges the alignment of the two blocks and returns a new alignment and consensus.
    - the first node is considered the anchor node, fixed in strandedness
    - the second node is appended to the first. The orientation and side depends on the edge.
    """
    new_node_ids = {k: n.id for k, n in new_nodes.items()}
    b1 = graph.blocks[edge.n1.bid]
    b2 = graph.blocks[edge.n2.bid]

    # make sure n1 is forward when merging
    if edge.n1.strand:
        b_left, b_right = b1, b2
    else:
        b_left, b_right = b2, b1

    # if they are not co-oriented then reverse complement the right block
    if edge.n1.strand != edge.n2.strand:
        b_right = b_right.reverse_complement()

    # assign the id equal to the id of the first block of the edge
    new_block_id = edge.n1.bid
    new_block = concatenate_alignments(
        b_left, b_right, node_map, new_node_ids, new_block_id
    )

    return new_block


def graph_merging_update_paths(
    graph: Pangraph, new_nodes: dict[int, Node], bid_left: int
):
    """
    Updates the paths in the graph: removes nodes corresponding to the two merged
    blocks and substitutes them with a new node for the merged block.
    """
    for pid, path in graph.paths.items():
        for i, nid in enumerate(path.nodes):
            if nid in new_nodes:
                old_node_bid = graph.nodes[nid].block_id
                if old_node_bid == bid_left:
                    path.nodes[i] = new_nodes[nid].id
                else:
                    path.nodes[i] = None
        # remove missing nodes
        path.nodes = [nid for nid in path.nodes if nid is not None]


def graph_merging_update_nodes(
    graph: Pangraph, new_nodes: dict[int, Node], bid_left: int
):
    """
    Updates the node dictionary in the graph: removes nodes corresponding
    to the two merged blocks and substitutes them with a new node for the merged block.
    """
    for nid, n in new_nodes.items():
        if graph.nodes[nid].block_id == bid_left:
            graph.nodes[n.id] = n
        del graph.nodes[nid]


def graph_merging_update(
    graph: Pangraph,
    new_block: Block,
    new_nodes: dict[int, Node],
    edge: Edge,
):
    """
    Updates the graph inplace after merging two blocks.
    """
    # delete old blocks and add new one
    del graph.blocks[edge.n1.bid]
    del graph.blocks[edge.n2.bid]
    graph.blocks[new_block.id] = new_block

    bid_left = edge.n1.bid

    # update paths
    graph_merging_update_paths(graph, new_nodes, bid_left)

    # delete old nodes and add new ones
    graph_merging_update_nodes(graph, new_nodes, bid_left)
