import utils as ut
from utils import Pangraph, Path, Block, Node, Substitution
from collections import Counter, defaultdict
from remove_nodes import remove_emtpy_nodes


# reconsensus function:
# - for blocks that:
#   - originate from a new merger
#   - have depth > 2
# - goes through each block and re-defines the consensus
#   - for mutations, changes any mutated position. Also update the alignment (this is straightforward)
#   - for any in/del present in > M/2 sites, appends it to the consensus
# - if the consensus has updated indels, then re-aligns all the sequences to the new consensus


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
    # reconstruct block sequences
    # NB: the "apply_edits_to_ref" function is the one defined in:
    # https://github.com/neherlab/pangraph/blob/d4272ece3ee326ca8228652667f9682f34d570a1/packages/pangraph/src/pangraph/edits.rs#L172
    seqs = {
        nid: ut.apply_edits_to_ref(edit, block.consensus)
        for nid, edit in block.alignment.items()
    }
    # update consensus
    block.consensus = consensus
    # re-align sequences
    # NB: map_variations is the function defined in:
    # https://github.com/neherlab/pangraph/blob/d4272ece3ee326ca8228652667f9682f34d570a1/packages/pangraph/src/align/map_variations.rs#L8
    # that uses nextalign to map the sequence to the consensus and extract the variations.
    # here I re-implemented it for simplicity using biopython.
    block.alignment = {
        nid: ut.map_variations(consensus, seq) for nid, seq in seqs.items()
    }
