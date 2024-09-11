import utils as ut
from utils import Pangraph, Path, Block, Node, Substitution
from collections import Counter, defaultdict

# remove empty nodes:
# - for blocks that originate from a new merger or split
# - goes through each block and removes any node with an empty sequence
# - also remove the node from the graph node list and paths

# reconsensus function:
# - for blocks that:
#   - originate from a new merger
#   - have depth > 2
# - goes through each block and re-defines the consensus
#   - for mutations, changes any mutated position
#   - for deletions, removes from consensus any
#   - for positions with many insertions, aligns them and make them surface on the consensus (MSA?)
#   - what if in/dels are nested?
# NB: only apply to **new mergers**! not to any updated blocks. Extract from merge promises.


def find_empty_nodes(graph: Pangraph, block_ids: list[int]) -> list[int]:
    """
    Finds nodes with empty sequences (full deletion) in the graph.
    It only checks specified blocks (the ones that were updated by a merger/split).
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
                # debug assert: single deletion
                assert len(edits.dels) == 1

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


def reconsensus_graph(graph: Pangraph, ids_updated_blocks: list[int]):
    """
    Applies the reconsensus operation to each updated block in the graph:
    - updates the block consensus following a merge
    - removes potentially empty nodes
    """
    node_ids_to_delete = []
    for block_id in ids_updated_blocks:
        block = graph.blocks[block_id]
        node_ids = reconsensus(block)
        node_ids_to_delete += node_ids

    # remove selected nodes from graph
    graph_remove_nodes(graph, node_ids_to_delete)


def reconsensus(block: Block) -> list[int]:
    """
    Performs the reconsensus operation inplace on a block.
    - if a position is mutated in > N/2 sites, it re-computes the consensus
    - if a position is deleted in > N/2 sites, it removes it from the consensus and saves changes as insertions
    - if a position has insertions in > N/2 sites, it computes a MSA for the insertions in the region, and potentially updates the consensus accordingly.
    """
    pass


def reconsensus_mutations(block: Block):
    """
    Re-computes the consensus for a block if a position is mutated in > N/2 sites,
    excluding deleted sites.
    """
    N = len(block.alignment)
    aln = block.alignment
    muts = defaultdict(Counter)
    # count mutations
    for nid, edit in aln.items():
        for s in edit.subs:
            muts[s.pos].update([s.alt])

    # change positions that are different in more than N>2 sites.
    changes = []
    for pos, ct in muts.items():
        alt, n = ct.most_common(1)[0]
        if n > N / 2:
            changes.append((pos, alt))

    # apply change
    for pos, alt in changes:
        # update consensus
        original = block.consensus[pos]
        bc = list(block.consensus)
        bc[pos] = alt
        block.consensus = "".join(bc)
        # change mutations
        for nid, edit in aln.items():
            subs = edit.subs
            subs_at_pos = list(filter(lambda s: s.pos == pos, subs))
            if len(subs_at_pos) == 0:
                subs.append(Substitution(pos, original))
                subs.sort(key=lambda s: s.pos)
            elif len(subs_at_pos) == 1:
                s = subs_at_pos[0]
                if s.alt == alt:
                    subs.remove(s)
            else:
                raise ValueError(
                    f"more than one substitution at site {pos}: {subs_at_pos}"
                )


def reconstruct_deletions(block: Block):
    """
    Removes a position from the consensus if it is deleted in > N/2 sites.
    """
    pass


def reconsensus_insertions(block: Block):
    """
    Computes a MSA for insertions in a region if they are present in > N/2 sites.
    """
    pass


def graph_remove_nodes(graph: Pangraph, node_ids: list[int]):
    """
    Deletes a list of nodes from the graph inplace. They are removed from the respective
    paths, from the respective block alignment dictionary, and from the node dictionary.
    """
