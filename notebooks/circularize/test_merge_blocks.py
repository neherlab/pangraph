from .utils import Block, Node, Path, Pangraph, Edit, Substitution, Insertion, Deletion
from .circularize_utils import SimpleNode, Edge
from .circularize import remove_transitive_edges
from .merge_blocks import (
    merge_blocks,
    find_node_pairings,
    concatenate_alignments,
    graph_merging_update_nodes,
    graph_merging_update_paths,
)
import pytest


@pytest.fixture
def block_1():
    #          0         1         2         3
    #          01234567890123456789012345678901
    # cons:    ACTATATTACGGCGATCGATCGATTACTCGCT
    #   n1:    ...G............................  l = 32
    #   n2:    .......|.....xxx................  l = 31
    #   n3:    ................................| l = 35
    return Block(
        id=1,
        consensus="ACTATATTACGGCGATCGATCGATTACTCGCT",
        alignment={
            1: Edit(ins=[], subs=[Substitution(3, "G")], dels=[]),
            2: Edit(ins=[Insertion(7, "AA")], subs=[], dels=[Deletion(13, 3)]),
            3: Edit(ins=[Insertion(32, "CCC")], subs=[], dels=[]),
        },
    )


@pytest.fixture
def block_2():
    #          0         1         2         3
    #          0123456789012345678901234567890
    # cons:    GATCTTAGGATCATCCCTATCATAGGAGTCG
    #   n4:    .........................xx....  l = 29
    #   n5:    ...T...........................  l = 31
    #   n6:   |xx.............................  l = 32
    return Block(
        id=2,
        consensus="GATCTTAGGATCATCCCTATCATAGGAGTCG",
        alignment={
            4: Edit(ins=[], subs=[], dels=[Deletion(25, 2)]),
            5: Edit(ins=[], subs=[Substitution(3, "T")], dels=[]),
            6: Edit(ins=[Insertion(0, "TTT")], subs=[], dels=[Deletion(0, 2)]),
        },
    )


@pytest.fixture
def block_3():
    #          0         1         2
    #          012345678901234567890
    # cons:    CTATTACTAGGGGGACCACTA
    #   n7:    ...............xx....  l = 19
    #   n8:    ...C.................  l = 21
    return Block(
        id=3,
        consensus="CTATTACTAGGGGGACCACTA",
        alignment={
            7: Edit(ins=[], subs=[], dels=[Deletion(15, 2)]),
            8: Edit(ins=[], subs=[Substitution(3, "C")], dels=[]),
        },
    )


@pytest.fixture
def graph_A(block_1, block_2, block_3):
    #      (0|32)      (32|61)     (61|0)
    # p1) (b1+|n1) -> (b2-|n4) -> (b3+|n7)  l=80
    #      (10|41)     (41|72)     (72|10)
    # p2) (b1+|n2) -> (b2-|n5) -> (b3+|n8)  l=83
    #      (5|40)      (40|5)
    # p3) (b2+|n6) -> (b1-|n3)              l=67

    paths = {
        1: Path(id=1, nodes=[1, 4, 7], L=80, circular=True),
        2: Path(id=2, nodes=[2, 5, 8], L=83, circular=True),
        3: Path(id=3, nodes=[6, 3], L=67, circular=True),
    }
    blocks = {b.id: b for b in [block_1, block_2, block_3]}
    nodes = {
        1: Node(id=1, block_id=1, path_id=1, position=(0, 32), strandedness=True),
        2: Node(id=2, block_id=1, path_id=2, position=(10, 41), strandedness=True),
        3: Node(id=3, block_id=1, path_id=3, position=(40, 5), strandedness=False),
        4: Node(id=4, block_id=2, path_id=1, position=(32, 61), strandedness=False),
        5: Node(id=5, block_id=2, path_id=2, position=(41, 72), strandedness=False),
        6: Node(id=6, block_id=2, path_id=3, position=(5, 40), strandedness=True),
        7: Node(id=7, block_id=3, path_id=1, position=(61, 0), strandedness=True),
        8: Node(id=8, block_id=3, path_id=2, position=(72, 10), strandedness=True),
    }
    return Pangraph(paths=paths, blocks=blocks, nodes=nodes)


def test_find_node_pairings_A(graph_A):
    n1 = SimpleNode(1, True)
    n2 = SimpleNode(2, False)
    edge = Edge(n1, n2)
    pairings, new_nodes = find_node_pairings(graph_A, edge)
    assert pairings == {1: 4, 2: 5, 3: 6, 4: 1, 5: 2, 6: 3}


@pytest.fixture
def block_1_revcomp():
    #          0         1         2         3
    #          01234567890123456789012345678901
    # cons:    AGCGAGTAATCGATCGATCGCCGTAATATAGT
    #   n1:    ............................C...  l = 32
    #   n2:    ................xxx......|......  l = 31
    #   n3:    |...............................  l = 35
    return Block(
        id=1,
        consensus="AGCGAGTAATCGATCGATCGCCGTAATATAGT",
        alignment={
            1: Edit(ins=[], subs=[Substitution(28, "C")], dels=[]),
            2: Edit(ins=[Insertion(25, "TT")], subs=[], dels=[Deletion(16, 3)]),
            3: Edit(ins=[Insertion(0, "GGG")], subs=[], dels=[]),
        },
    )


@pytest.fixture
def block_2_revcomp():
    #          0         1         2         3
    #          0123456789012345678901234567890
    # cons:    CGACTCCTATGATAGGGATGATCCTAAGATC
    #   n4:    ....xx.........................  l = 29
    #   n5:    ...........................A...  l = 31
    #   n6:    .............................xx| l = 32
    return Block(
        id=2,
        consensus="CGACTCCTATGATAGGGATGATCCTAAGATC",
        alignment={
            4: Edit(ins=[], subs=[], dels=[Deletion(4, 2)]),
            5: Edit(ins=[], subs=[Substitution(27, "A")], dels=[]),
            6: Edit(ins=[Insertion(31, "AAA")], subs=[], dels=[Deletion(29, 2)]),
        },
    )


def test_reverse_complement_1(block_1, block_1_revcomp):
    rev_block = block_1.reverse_complement()
    assert rev_block == block_1_revcomp


def test_reverse_complement_2(block_2, block_2_revcomp):
    rev_block = block_2.reverse_complement()
    assert rev_block == block_2_revcomp


@pytest.fixture
def expected_concat():
    #          0         1         2         3         4         5         6
    #          012345678901234567890123456789012345678901234567890123456789012
    # cons:    ACTATATTACGGCGATCGATCGATTACTCGCTCGACTCCTATGATAGGGATGATCCTAAGATC
    #   n1:    ...G................................xx.........................  l = 32 + 29
    #   n2:    .......|.....xxx...........................................A...  l = 31 + 31
    #   n3:    ................................|............................xx| l = 35 + 32
    return Block(
        id=1,
        consensus="ACTATATTACGGCGATCGATCGATTACTCGCTCGACTCCTATGATAGGGATGATCCTAAGATC",
        alignment={
            -9124324634939260889: Edit(
                ins=[],
                subs=[Substitution(3, "G")],
                dels=[Deletion(36, 2)],
            ),
            -8785834526545769370: Edit(
                ins=[Insertion(7, "AA")],
                subs=[Substitution(59, "A")],
                dels=[Deletion(13, 3)],
            ),
            893596598396622086: Edit(
                ins=[Insertion(32, "CCC"), Insertion(63, "AAA")],
                subs=[],
                dels=[Deletion(61, 2)],
            ),
        },
    )


@pytest.fixture
def new_node_ids():
    return {
        1: -9124324634939260889,
        2: -8785834526545769370,
        3: 893596598396622086,
        4: -9124324634939260889,
        5: -8785834526545769370,
        6: 893596598396622086,
    }


def test_concatenate_blocks(graph_A, block_1, block_2, expected_concat, new_node_ids):
    n1 = SimpleNode(1, True)
    n2 = SimpleNode(2, False)
    edge = Edge(n1, n2)
    pairings, new_nodes = find_node_pairings(graph_A, edge)
    new_block_id = 1
    block = concatenate_alignments(
        block_1, block_2.reverse_complement(), pairings, new_node_ids, new_block_id
    )
    assert block == expected_concat


@pytest.fixture
def expected_graph(expected_concat, new_node_ids):
    #       (0|-----------|61)     (61|0)
    # p1) (b1+|-----------|n1) -> (b3+|n7)  l=80
    #      (10|-----------|72)     (72|10)
    # p2) (b1+|-----------|n2) -> (b3+|n8)  l=83
    #       (5|-----------|5)
    # p3) (b1-|-----------|n3)              l=67
    blocks = {
        1: expected_concat,
        3: Block(
            id=3,
            consensus="CTATTACTAGGGGGACCACTA",
            alignment={
                7: Edit(ins=[], subs=[], dels=[Deletion(15, 2)]),
                8: Edit(ins=[], subs=[Substitution(3, "C")], dels=[]),
            },
        ),
    }

    nodes = {
        new_node_ids[1]: Node(
            id=new_node_ids[1],
            block_id=1,
            path_id=1,
            position=(0, 61),
            strandedness=True,
        ),
        new_node_ids[2]: Node(
            id=new_node_ids[2],
            block_id=1,
            path_id=2,
            position=(10, 72),
            strandedness=True,
        ),
        new_node_ids[3]: Node(
            id=new_node_ids[3],
            block_id=1,
            path_id=3,
            position=(5, 5),
            strandedness=False,
        ),
        7: Node(id=7, block_id=3, path_id=1, position=(61, 0), strandedness=True),
        8: Node(id=8, block_id=3, path_id=2, position=(72, 10), strandedness=True),
    }

    paths = {
        1: Path(id=1, nodes=[new_node_ids[1], 7], L=80, circular=True),
        2: Path(id=2, nodes=[new_node_ids[2], 8], L=83, circular=True),
        3: Path(id=3, nodes=[new_node_ids[3]], L=67, circular=True),
    }
    print(nodes)
    return Pangraph(paths=paths, blocks=blocks, nodes=nodes)


def test_graph_merging_update_paths(graph_A, expected_graph, new_node_ids):
    n1 = SimpleNode(1, True)
    n2 = SimpleNode(2, False)
    edge = Edge(n1, n2)
    pairings, new_nodes = find_node_pairings(graph_A, edge)
    graph_merging_update_paths(graph_A, new_nodes, bid_left=1)
    assert graph_A.paths == expected_graph.paths


def test_graph_merging_update_nodes(graph_A, expected_graph, new_node_ids):
    n1 = SimpleNode(1, True)
    n2 = SimpleNode(2, False)
    edge = Edge(n1, n2)
    pairings, new_nodes = find_node_pairings(graph_A, edge)
    graph_merging_update_nodes(graph_A, new_nodes, bid_left=1)
    assert graph_A.nodes == expected_graph.nodes


def test_merge_blocks(graph_A, expected_graph):
    n1 = SimpleNode(1, True)
    n2 = SimpleNode(2, False)
    edge = Edge(n1, n2)
    merge_blocks(graph_A, edge)
    assert graph_A.nodes == expected_graph.nodes
    assert graph_A.blocks == expected_graph.blocks
    assert graph_A.paths == expected_graph.paths


def test_remove_transitive_edges(graph_A, expected_graph):
    remove_transitive_edges(graph_A)
    assert graph_A.nodes == expected_graph.nodes
    assert graph_A.blocks == expected_graph.blocks
    assert graph_A.paths == expected_graph.paths
