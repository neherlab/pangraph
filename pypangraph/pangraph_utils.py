from collections import defaultdict, Counter
import networkx as nx
import numpy as np


def pangraph_to_networkx(pan, only_blocks=None):
    """Given a pangraph object (and optionally a list of block ids) this function
    returns a graph in which nodes are blocks, and edges connect two blocks
    that occurr one after the other.

    The optional argument `only_blocks` must be a list of block ids. If passed
    then the graph is built by excluding from the paths any block that is not
    part of this list.

    The returned graph object has one additional property: paths.
    This is a dictionary {strain -> list of block ids} and blocks is a list
    of all block ids present in the graph
    """

    # express the pangraph as a dictionary of strains -> list of block ids
    paths_dict = pan.to_paths_dict()

    # create graph blocks and paths
    G_blocks = pan.block_ids()
    G_paths = paths_dict

    # if a block restriction is requested, then restrict the objects
    if only_blocks is not None:
        G_blocks = only_blocks
        for strain in G_paths:
            G_paths[strain] = [bl for bl in G_paths[strain] if bl in G_blocks]

    # initialize empty graph
    G = nx.Graph()

    # cycle through all paths
    for strain, path in G_paths.items():

        # add all blocks from the path as nodes
        G.add_nodes_from(path)

        # cycle through consecutive block pairs
        for x, y in zip(path, np.roll(path, -1)):
            # add the edge to the graph
            G.add_edge(x, y)

    # append extra property
    G.paths = G_paths

    # return graph
    return G


def node_attribute(G, pan, attr):
    """Evaluates a specific attribute of nodes in the graph, according to the
    `attr` argument. Possible values are:
    - n_strains : the number of strains in which a given block occurrs at least once
    - n_occurrences : the total number of occurrences of the block, possibly with repetitions.
    - seq_len : the sequence length of the block
    - n_mut : the average number of mutations per block occurrence identified by pangraph

    Args:
    - G (networkx Graph) : this must be the graph created by the pangraph 'to_networkx'
        function.
    - pan (Pangraph) : the orginal pangraph object with which the graph was created
    - attr (string) : node attribute to evaluate. Must have one of the values specified above.
    """

    # number of different strains that possess a block
    if attr == "n_strains":
        strain_counter = Counter()
        for path in G.paths.values():
            # count every unique block once.
            strain_counter.update(np.unique(path))
        # dictionary {block_id -> n. strains}
        return dict(strain_counter)

    # number of total repetitions of a block
    if attr == "n_occurrences":
        block_counter = Counter()
        for path in G.paths.values():
            # count every single block occurrence
            block_counter.update(path)
        # dictionary {block_id -> n. block occurrences}
        return dict(block_counter)

    # length of block in nucleotides
    if attr == "seq_len":
        length_dict = {bl: len(pan.blocks[bl]) for bl in G.nodes()}
        # dictionary {block_id -> len block in nucleotides}
        return length_dict

    # n. SNPs in block
    if attr == "n_mut":
        n_mut = {}
        for bl in G.nodes():
            al_diff = pan.blocks[bl].align_diff
            n_mut_bl = [len(item["mutate"]) for item in al_diff.values()]
            n_mut[bl] = np.mean(n_mut_bl)
        return n_mut

    raise ValueError(
        "the attribute must be one of the following strings:\n \
     n_strains , n_occurrences , seq_len , n_mut"
    )


def edge_attribute(G, pan, attr):
    """Evaluates a specific attribute of edges in the graph, according to the
    `attr` argument. Possible values are:
    - n_occurrences : the total number of occurrences of the edge.
    - n_hidden_blocks : the total number of hidden blocks in between two edges.
    - n_hidden_blocks_occurrences : the total number of block occurrences that
        have been excluded, that would otherwise lay between two blocks that remain.
        Note that this will also count block repetitions.

    Nb: it returns a dictionary in the form {(x,y) -> value}, where for ease of
    search both orders (x,y) and (y,x) are present. This is for compatibility
    with networkx, in which node order in edges is not alphabetic.

    Args:
    - G (networkx Graph) : this must be the graph created by the pangraph 'to_networkx'
        function.
    - pan (Pangraph) : the orginal pangraph object with which the graph was created
    - attr (string) : edge attribute to evaluate. Must have one of the values specified above.
    """

    # n. of times the two blocks occurr one after the other
    if attr == "n_occurrences":
        occ_dict = defaultdict(int)
        for strain, path in G.paths.items():
            # for every path
            for x, y in zip(path, np.roll(path, -1)):
                # count every successive block occurrence once
                occ_dict[(x, y)] += 1
                occ_dict[(y, x)] += 1
        return dict(occ_dict)

    # n. of hidden blocks that have been compressed on the edge.
    # Nb: if the same block is repeated, it is counted multiple times.
    if attr == "n_hidden_blocks_occurrences":

        path_dict = pan.to_paths_dict()
        kept_blocks = G.nodes()  # list of blocks that are not hidden

        # counter for block occurrences on an edge
        edge_compression = defaultdict(int)

        for path in path_dict.values():
            # names of blocks that are kept, and number of hidden blocks in between
            ids_kept, n_excluded = [], [0]

            for bl in path:
                # if not hidden, then go to the next edge
                if bl in kept_blocks:
                    ids_kept.append(bl)
                    n_excluded.append(0)
                else:
                    # if hidden, add one more hidden block
                    n_excluded[-1] += 1

            # compensate for periodic boundary conditions
            n_excluded[0] += n_excluded[-1]
            n_excluded = n_excluded[:-1]

            # turn the results in a dictionary whose keys are pairs of nodes.
            # Nb: every count is repeated twice because nodes are ordered in both ways.
            for n_edge, (x, y) in enumerate(zip(np.roll(ids_kept, 1), ids_kept)):
                edge_compression[(x, y)] += n_excluded[n_edge]
                edge_compression[(y, x)] += n_excluded[n_edge]

        return dict(edge_compression)

    if attr == "n_hidden_blocks":

        path_dict = pan.to_paths_dict()
        kept_blocks = G.nodes()  # list of blocks that are not hidden

        # counter for block occurrences on an edge
        edge_hiddenblocks = defaultdict(list)

        for path in path_dict.values():
            # names of blocks that are kept, and names of hidden blocks in between
            ids_kept, n_excluded = [], [[]]

            for bl in path:
                # if not hidden, then go to the next edge
                if bl in kept_blocks:
                    ids_kept.append(bl)
                    n_excluded.append([])
                else:
                    # if hidden, add one more hidden block
                    n_excluded[-1].append(bl)

            # compensate for periodic boundary conditions
            n_excluded[0] += n_excluded[-1]
            n_excluded = n_excluded[:-1]

            # turn the results in a dictionary whose keys are pairs of nodes.
            # Nb: every count is repeated twice because nodes are ordered in both ways.
            for n_edge, (x, y) in enumerate(zip(np.roll(ids_kept, 1), ids_kept)):
                edge_hiddenblocks[(x, y)] += n_excluded[n_edge]
                edge_hiddenblocks[(y, x)] += n_excluded[n_edge]

        # turn from list of nodes to count of unique values
        edge_hiddenblocks = {
            k: len(np.unique(edge_hiddenblocks[k])) for k in edge_hiddenblocks
        }

        return edge_hiddenblocks

    raise ValueError(
        "the attribute must be one of the following strings:\n \
     n_occurrences"
    )
