# here is an overview of how I imagine the algorithm will look like,
# for the graph merging part. This is the function that is called
# when a node of the guide tree is visited.
#
#   left graph        right graph
#             \      /
#              \    /
#           merged graph
#
# In the visited node we will call the merge_graphs function, which
# will take the two graphs and merge them into a single one. The merged
# graph is passed to the parent node.

from dataclasses import dataclass
from pangraph_classes import Block, Path, Pangraph


class Alignment:
    # this is the alignment class that you already implemented in rust,
    # to deal with paf files. It contains information on the homology
    # between two blocks.
    pass


def merge_graphs(left_graph: Pangraph, right_graph: Pangraph) -> Pangraph:

    # put the two graphs in a single one, by simply joining
    # the two sets of blocks and paths. No merging is pefrormed
    graph = graph_join(left_graph, right_graph)

    # iteratively try to merge homologous regions in blocks.
    # We map the consensus sequences of blocks to each other and merge
    # blocks in which we find matches. We iterate this until no more merging
    # is possible.
    while True:
        graph, has_changed = self_merge(graph)
        # stop when no more mergers are possible
        if not has_changed:
            break

    return graph


def graph_join(left_graph: Pangraph, right_graph: Pangraph) -> Pangraph:
    # simply join the two sets of blocks and paths
    join_blocks = left_graph.blocks | right_graph.blocks
    join_paths = left_graph.paths | right_graph.paths
    graph = Pangraph(blocks=join_blocks, paths=join_paths)
    return graph


def self_merge(graph: Pangraph) -> Pangraph:

    # use minimap2 or other aligners to find matches between the consensus
    # sequences of the blocks
    matches = find_matches(graph.blocks)

    # split matches:
    # - whenever an alignment contains an in/del longer than the threshold length
    #   (parameter - default 100 bp) we want to split the alignment in two)
    matches = split_matches(matches)

    # filter matches:
    # - calculate energy and keep only matches with E < 0
    # - sort them by energy
    # - discard incompatible matches (the ones that have overlapping regions)
    matches = filter_matches(matches)

    # complex function: takes the list of desired matches and the two
    # graphs. It splits blocks and re-weaves paths through them. Paths
    # will not be updated after this.
    # The function does not perform any merging, but returns a list of
    # blocks that should be merged.
    # Nb: this function should already take care of:
    # - creating new nodes
    # - substituting the old nodes with the new ones in paths
    # - splitting blocks and updating the node ids.
    # - adding the blocks that do not need merging to the preliminary graph
    # - return the set of blocks that should be merged
    graph, mergers = reweave_graph(graph, matches)

    # this can be parallelized
    merged_blocks = {}
    for merger in mergers:
        # merge two blocks together and return a
        block_id = merger.merged_id
        merged_block = perform_merger(merger)
        # append to the block dictionary
        merged_blocks[block_id] = merged_blocks

    # add the new blocks to the graph
    graph.blocks |= merged_blocks

    # we might need this step for some final updates and consistency checks.
    # TODO: here we could also take care of transitive edges, which is useful
    # in the case of circular paths.
    graph = consolidate(graph)

    return graph


def find_matches(blocks: dict[str, Block]) -> list[Alignment]:
    """This function calls an aligner (default: minimap2) to find matches
    between the consensus sequences of the blocks. It returns a list of
    alignment objects."""

    seqs = [block.consensus for block in blocks]
    alignments = aligner(seqs)
    return alignments


def aligner(seqs: list[str]) -> list[Alignment]:
    # use minimap2 or other aligners to find matches between the consensus
    # sequences of the blocks.
    pass


def split_matches(matches: list[Alignment]) -> list[Alignment]:
    # TODO: split the alignments whenever an alignment contains an in/del
    # longer than the threshold length (parameter - default 100 bp).
    # There are some edge cases that we might need to consider here
    # (eg. 1000M 200 D 1M 100I 200 M)
    # returns the updated list of alignments.
    pass


def filter_matches(matches: list[Alignment]) -> list[Alignment]:
    # - evaluates the energy of the alignments
    # - keeps only matches with E < 0
    # - sorts them by energy
    # - discards incompatible matches (the ones that have overlapping regions)
    energies = [alignment_energy(match) for match in matches]

    # keep only matches with E < 0
    matches = [match for match, energy in zip(matches, energies) if energy < 0]

    # sort them by increasing value of energy
    matches = [match for energy, match in sorted(zip(energies, matches))]

    # discard incompatible matches
    keep_matches = []
    for match in matches:
        if compatible_match(match, keep_matches):
            keep_matches.append(match)

    return keep_matches


def alignment_energy(alignment: Alignment) -> float:
    # TODO: calculate the energy of the alignment
    # this is a function of alignment length, identity, and whether it will
    # create a block split.
    pass


def compatible_match(match: Alignment, matches: list[Alignment]) -> bool:
    # TODO: check that the new match does not overlap with any of the
    # matches that we already have.
    pass


@dataclass()
class Merger:
    merged_block_id: str
    anchor_block: Block
    shallow_block: Block
    alignment: Alignment


def reweave_graph(
    graph: Pangraph, matches: list[Alignment]
) -> tuple[Pangraph, list[Merger]]:
    # TODO: complex function. I will expand more on this, but it should:
    # - create a new graph with a copy of the paths, and only the blocks
    #   that do no undergo any merging.
    # - for each of the suggested matches:
    #   - split the block and create new nodes (including ones for mergers)
    #   - add the blocks that should not be processed further to the new graph
    #   - update the path with the new nodes
    #   - for blocks that should be merged, create a merger object and add it
    #     to the list of mergers.
    # TODO: this function should also update the position of all of the nodes!
    pass


def perform_merger(merger: Merger) -> Block:
    # TODO: this part is basically contained in the block_operations notebook.
    # in the merger we add all sequences from the shallow block to the deep
    # block, and return the updated block.

    # TODO: potentially here later we can also take care of consolidating the
    # consensus sequence, updating it with indels or substitutions that might
    # have become present in the majority of sequences.
    pass


def consolidate(graph: Pangraph) -> Pangraph:
    # TODO: final updates and consistency checks.
    # - we can take care of removing _transitive edges_, i.e.
    #   paris of blocks that are always adjacent and connected in the same
    #   way. We have an algorithm to quickly check for this. This is useful when
    #   dealing with circular paths.
    # - we can also add optional consistency checks, to make sure that sequence
    #   lengths and path lengths are conserved, and optionally check that we can.
    #   reconstruct the full genomes exactly.
    pass
