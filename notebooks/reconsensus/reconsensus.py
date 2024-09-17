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
#   - for mutations, changes any mutated position. Also update the alignment (this is straightforward)
#   - for any in/del present in > M/2 sites, appends it to the consensus
# - if the consensus has updated indels, then re-aligns all the sequences to the new consensus


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


def reconsensus_graph(graph: Pangraph, ids_updated_blocks: list[int]):
    """
    Applies the reconsensus operation to each updated block in the graph:
    - updates the block consensus following a merge
    - removes potentially empty nodes
    """
    # remove selected nodes from graph
    remove_emtpy_nodes(graph, ids_updated_blocks)

    for block_id in ids_updated_blocks:
        block = graph.blocks[block_id]
        reconsensus(block)


def reconsensus(block: Block):
    """
    Performs the reconsensus operation inplace on a block.
    - if a position is mutated in > N/2 sites, it adds the mutation to the consensus and updates the alignment.
    - if an in/del is present in > N/2 sites, it adds it to the consensus and re-aligns the sequences to the updated consensus.
    """
    reconsensus_mutations(block)
    ins = majority_insertions(block)
    dels = majority_deletions(block)
    if len(ins) > 0 or len(dels) > 0:
        cons = block.consensus
        cons = apply_indels(cons, dels, ins)
        update_block_consensus(block, cons)


def reconsensus_mutations(block: Block):
    """
    Re-computes the consensus for a block if a position is mutated in > N/2 sites.
    """
    N = block.depth()
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


def majority_deletions(block: Block) -> list[int]:
    """
    Returns a list of positions to be removed from the consensus, because they are deleted in > N/2 sites.
    """
    N = block.depth()
    n_dels = Counter()
    for nid, edit in block.alignment.items():
        for d in edit.dels:
            # for each deleted position, increment the counter
            n_dels.update(range(d.pos, d.pos + d.length))
    # return the positions that are deleted in more than N/2 sites
    return [pos for pos, n in n_dels.items() if n > N / 2]


def majority_insertions(block: Block) -> list[tuple[int, str]]:
    """
    Returns a list of insertions to be added to the consensus, because they are inserted in > N/2 sites.
    """
    N = block.depth()
    n_ins = Counter()
    for nid, edit in block.alignment.items():
        n_ins.update([(i.pos, i.ins) for i in edit.ins])
    # return the positions that are inserted in more than N/2 sites
    return [(pos, ins) for (pos, ins), n in n_ins.items() if n > N / 2]


def apply_indels(seq: str, dels: list[int], ins: list[tuple[int, str]]) -> str:
    """
    Updates the consensus sequence with the deletions and insertions.
    """
    seq = list(seq)
    for pos in dels:
        seq[pos] = ""
    for pos, ins in sorted(ins, key=lambda x: x[0], reverse=True):
        seq.insert(pos, ins)
    return "".join(seq)


def update_block_consensus(block: Block, consensus: str):
    """
    updates the consensus sequence of the block and re-aligns the sequences to the new consensus.
    """
    pass
