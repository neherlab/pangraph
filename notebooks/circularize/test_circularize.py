import pytest

from .utils import Pangraph, Path, Block, Node, Edit, Substitution, Insertion, Deletion
from .circularize import (
    Edge,
    SimpleNode,
    block_depths,
    edges_count,
    transitive_edges,
)


def test_SimpleNode():
    n1 = SimpleNode(1, True)
    n2 = SimpleNode(1, False)
    n3 = SimpleNode(2, True)
    n4 = SimpleNode(2, False)
    assert n1 != n2
    assert n1 == ~n2
    assert ~n1 == n2
    assert ~n1 != n3
    assert n1 != n3
    assert n1 != ~n3
    assert n4 != n3
    assert n4 == ~n3


def test_Edge():
    n1 = SimpleNode(1, True)
    n2 = SimpleNode(2, False)

    e1 = Edge(n1, n2)
    e2 = Edge(n2, n1)
    e3 = Edge(~n2, ~n1)

    assert e1 != e2
    assert e1 == e3

    x = {e1: 1, e2: 2}
    assert e3 in x

    assert set([e1, e2, e3]) == {e1, e2}


@pytest.fixture
def input_graph():
    # create the following graph with circular paths:
    # a) 1+(10) -> 2+(20) -> 3+(30) -> 4+(40)
    # b) 1+(11) -> 2-(21) -> 2+(22) -> 3+(31) -> 4+(41)
    # c) 1+(12) -> 2+(23) -> 3-(32) -> 4+ (42)
    # d) 1+(13) -> 3-(33) -> 2+(24) -> 3-(34) -> 4+(43)
    # f) 4- (44) -> 3-(35) -> 2-(25) -> 1-(14)
    paths = {
        0: Path(0, nodes=[10, 20, 30, 40], L=None, circular=True),
        1: Path(1, nodes=[11, 21, 22, 31, 41], L=None, circular=True),
        2: Path(2, nodes=[12, 23, 32, 42], L=None, circular=True),
        3: Path(3, nodes=[13, 33, 24, 34, 43], L=None, circular=True),
        4: Path(4, nodes=[44, 35, 25, 14], L=None, circular=True),
    }
    nodes = {
        10: Node(10, 1, 0, (), True),
        20: Node(20, 2, 0, (), True),
        30: Node(30, 3, 0, (), True),
        40: Node(40, 4, 0, (), True),
        11: Node(11, 1, 1, (), True),
        21: Node(21, 2, 1, (), False),
        22: Node(22, 2, 1, (), True),
        31: Node(31, 3, 1, (), True),
        41: Node(41, 4, 1, (), True),
        12: Node(12, 1, 2, (), True),
        23: Node(23, 2, 2, (), True),
        32: Node(32, 3, 2, (), False),
        42: Node(42, 4, 2, (), True),
        13: Node(13, 1, 3, (), True),
        33: Node(33, 3, 3, (), False),
        24: Node(24, 2, 3, (), True),
        34: Node(34, 3, 3, (), False),
        43: Node(43, 4, 3, (), True),
        44: Node(44, 4, 4, (), False),
        35: Node(35, 3, 4, (), False),
        25: Node(25, 2, 4, (), False),
        14: Node(14, 1, 4, (), False),
    }

    def mock_aln(nids):
        return {nid: Edit(None, None, None) for nid in nids}

    blocks = {
        1: Block(
            id=1,
            consensus=None,
            alignment=mock_aln([10, 11, 12, 13, 14]),
        ),
        2: Block(
            id=2,
            consensus=None,
            alignment=mock_aln([20, 21, 22, 23, 24, 25]),
        ),
        3: Block(
            id=3,
            consensus=None,
            alignment=mock_aln([30, 31, 32, 33, 34, 35]),
        ),
        4: Block(
            id=4,
            consensus=None,
            alignment=mock_aln([40, 41, 42, 43, 44]),
        ),
    }
    return Pangraph(paths=paths, nodes=nodes, blocks=blocks)


def test_block_depths(input_graph):
    bc = block_depths(input_graph)
    assert bc == {1: 5, 2: 6, 3: 6, 4: 5}


def test_edges_count(input_graph):
    ec = edges_count(input_graph)
    n1 = SimpleNode(1, True)
    n2 = SimpleNode(2, True)
    n3 = SimpleNode(3, True)
    n4 = SimpleNode(4, True)
    assert ec[Edge(n1, n2)] == 3
    assert ec[Edge(n1, ~n2)] == 1
    assert ec[Edge(n2, n3)] == 3
    assert ec[Edge(n2, ~n3)] == 2
    assert ec[Edge(~n2, n2)] == 1
    assert Edge(n2, ~n2) not in ec
    assert ec[Edge(n3, n4)] == 3
    assert ec[Edge(~n3, n4)] == 2
    assert ec[Edge(n4, n1)] == 5


def test_transitive_edges_A(input_graph):
    te = transitive_edges(input_graph)
    n1 = SimpleNode(1, True)
    n4 = SimpleNode(4, True)
    assert te == [Edge(n4, n1)]


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
            3: Edit(ins=[Insertion(31, "CCC")], subs=[], dels=[]),
        },
    )


@pytest.fixture
def graph_B(block_1):
    # single block: no mergings
    #      (0|32)
    # p1) (b1+|n1)  l=32
    #      (0|31)
    # p2) (b1+|n2)  l=31
    #      (0|35)
    # p3) (b1-|n3)  l=35
    paths = {
        1: Path(id=1, nodes=[1], L=32, circular=True),
        2: Path(id=2, nodes=[2], L=31, circular=True),
        3: Path(id=3, nodes=[3], L=35, circular=True),
    }
    blocks = {1: block_1}
    nodes = {
        1: Node(id=1, block_id=1, path_id=1, position=(0, 32), strandedness=True),
        2: Node(id=2, block_id=1, path_id=2, position=(0, 31), strandedness=True),
        3: Node(id=3, block_id=1, path_id=3, position=(0, 35), strandedness=False),
    }
    return Pangraph(paths=paths, blocks=blocks, nodes=nodes)


def test_transitive_edges_B(graph_B):
    te = transitive_edges(graph_B)
    assert te == []
