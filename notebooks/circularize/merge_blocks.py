from utils import Pangraph, Block, reverse_complement
from circularize_utils import Edge, SimpleNode


def merge_blocks(graph: Pangraph, edge: Edge):
    """
    Merges the two blocks from the same edge into a single block,
    modifying the graph in place.
    """

    # first is anchor block, the second is the block to append
    # the anchor block gets extended, on left or right.
    # while the reference block gets deleted.
    # TODO: how to update node id with hash?

    edge = orient_merging_edge(graph, edge)
    # 1) find node pairs
    node_map = find_node_pairings(graph, edge)
    # 2) merge alignment
    pass


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
    D1, D2 = b1.depth(), b2.depth()
    assert D1 == D2, "merging requires blocks to have the same depth"
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
    Given a specific edge between blocks, returns pairings for all nodes in the blocks.
    """
    node_map = {}
    for path_id, path in graph.paths.items():
        N = len(path.nodes)  # number of nodes
        I = N if path.circular else N - 1  # number of edges
        for i in range(I):
            nid1 = path.nodes[i]
            nid2 = path.nodes[(i + 1) % N]
            n1 = graph.nodes[nid1]
            n2 = graph.nodes[nid2]
            bid1, bid2 = n1.block_id, n2.block_id
            s1, s2 = n1.strandedness, n2.strandedness
            sn1 = SimpleNode(bid1, s1)
            sn2 = SimpleNode(bid2, s2)
            e = Edge(sn1, sn2)
            if edge == e:
                node_map[nid1] = nid2
                node_map[nid2] = nid1
    return node_map


def concatenate_alignments(bl1, bl2, node_map):
    """Concatenates two blocks, with block id and node ids of the first block."""
    seq1 = bl1.consensus
    seq2 = bl2.consensus
    aln1 = bl1.alignment
    aln2 = bl2.alignment
    L1 = bl1.consensus_len()

    seq = seq1 + seq2
    aln = {}
    for nid1, e1 in aln1.items():
        nid2 = node_map[nid1]
        e2 = aln2[nid2]
        aln[nid1] = e1.concat(e2.shift(L1))
    return Block(bl1.id, seq, aln)


def merge_alignment(graph: Pangraph, edge: Edge, node_map: dict[int, int]):
    """
    Merges the alignment of the two blocks and returns a new alignment and consensus.
    - the first node is considered the anchor node, fixed in strandedness
    - the second node is appended to the first. The orientation and side depends on the edge.
    """
    pass
