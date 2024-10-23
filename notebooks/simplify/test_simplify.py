import pytest

from .utils import Pangraph, Path, Block, Node, Edit, Substitution, Insertion, Deletion
from .simplify import simplify


@pytest.fixture
def block_A():
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
def block_B():
    #          0         1         2         3
    #          01234567890123456789012345678901
    # cons:    CATGCTACGCTACGCATTATCGATCGCATCGA
    #   n4:    ..........G.....................  l = 32
    #   n5:    .............xxx................  l = 29
    #   n6:    ................................|  l = 35
    return Block(
        id=1,
        consensus="CATGCTACGCTACGCATTATCGATCGCATCGA",
        alignment={
            4: Edit(ins=[], subs=[Substitution(10, "G")], dels=[]),
            5: Edit(ins=[], subs=[], dels=[Deletion(13, 3)]),
            6: Edit(ins=[Insertion(32, "AAA")], subs=[], dels=[]),
        },
    )


@pytest.fixture
def block_C():
    #          0         1
    #          01234567890123456
    # cons:    ACGTGTACTAGTACTGC
    #   n7:    ................. l = 17
    #   n8:    ............C.... l = 17
    return Block(
        id=1,
        consensus="ACGTGTACTAGTACTGC",
        alignment={
            7: Edit(ins=[], subs=[], dels=[]),
            8: Edit(ins=[], subs=[Substitution(12, "C")], dels=[]),
        },
    )


nid11 = -956790416118992873
nid12 = -550362699546468170
new_block_id = 2


@pytest.fixture
def block_AB():
    #          0         1         2         3         4         5         6
    #          01234567890123456789012345678901 23456789012345678901234567890123
    # cons:    ACTATATTACGGCGATCGATCGATTACTCGCT CATGCTACGCTACGCATTATCGATCGCATCGA
    #  n10:    ...G............................ ..........G.....................
    #  n11:    .......|.....xxx................ .............xxx................

    return Block(
        id=new_block_id,
        consensus="ACTATATTACGGCGATCGATCGATTACTCGCTCATGCTACGCTACGCATTATCGATCGCATCGA",
        alignment={
            nid11: Edit(
                ins=[], subs=[Substitution(3, "G"), Substitution(42, "G")], dels=[]
            ),
            nid12: Edit(
                ins=[Insertion(7, "AA")],
                subs=[],
                dels=[Deletion(13, 3), Deletion(45, 3)],
            ),
        },
    )


@pytest.fixture
def nodes():
    return {
        1: Node(id=1, block_id=1, path_id=1, strandedness=True, position=(0, 32)),
        2: Node(id=2, block_id=1, path_id=2, strandedness=True, position=(0, 31)),
        3: Node(id=3, block_id=1, path_id=3, strandedness=True, position=(0, 35)),
        4: Node(id=4, block_id=2, path_id=1, strandedness=True, position=(32, 64)),
        5: Node(id=5, block_id=2, path_id=2, strandedness=True, position=(31, 60)),
        6: Node(id=6, block_id=2, path_id=3, strandedness=True, position=(35, 0)),
        7: Node(id=7, block_id=3, path_id=1, strandedness=True, position=(64, 0)),
        8: Node(id=8, block_id=3, path_id=2, strandedness=False, position=(60, 0)),
    }


@pytest.fixture
def graph(nodes, block_A, block_B, block_C):
    # n1+ -> n4+ -> n7+
    # n2+ -> n5+ -> n8-
    # n3+ -> n6-
    return Pangraph(
        nodes=nodes,
        blocks={1: block_A, 2: block_B, 3: block_C},
        paths={
            1: Path(id=1, name="pathA", nodes=[1, 4, 7], L=81, circular=True),
            2: Path(id=2, name="pathB", nodes=[2, 5, 8], L=77, circular=True),
            3: Path(id=3, name="pathC", nodes=[3, 6], L=70, circular=True),
        },
    )


@pytest.fixture
def expected_graph(block_AB, block_C):
    # n1+|n4+ -> n7+
    # n2+|n5+ -> n8-
    return Pangraph(
        nodes={
            nid11: Node(
                id=nid11,
                block_id=new_block_id,
                path_id=1,
                strandedness=True,
                position=(0, 64),
            ),
            nid12: Node(
                id=nid12,
                block_id=new_block_id,
                path_id=2,
                strandedness=True,
                position=(0, 60),
            ),
            7: Node(id=7, block_id=3, path_id=1, strandedness=True, position=(64, 0)),
            8: Node(id=8, block_id=3, path_id=2, strandedness=False, position=(60, 0)),
        },
        blocks={new_block_id: block_AB, 3: block_C},
        paths={
            1: Path(id=1, name="pathA", nodes=[nid11, 7], L=81, circular=True),
            2: Path(id=2, name="pathB", nodes=[nid12, 8], L=77, circular=True),
        },
    )


def test_remove_path(graph):
    graph.remove_path(1)

    assert graph.paths == {
        2: Path(id=2, name="pathB", nodes=[2, 5, 8], L=77, circular=True),
        3: Path(id=3, name="pathC", nodes=[3, 6], L=70, circular=True),
    }
    assert graph.nodes == {
        2: Node(id=2, block_id=1, path_id=2, strandedness=True, position=(0, 31)),
        3: Node(id=3, block_id=1, path_id=3, strandedness=True, position=(0, 35)),
        5: Node(id=5, block_id=2, path_id=2, strandedness=True, position=(31, 60)),
        6: Node(id=6, block_id=2, path_id=3, strandedness=True, position=(35, 0)),
        8: Node(id=8, block_id=3, path_id=2, strandedness=False, position=(60, 0)),
    }
    assert graph.blocks == {
        1: Block(
            id=1,
            consensus="ACTATATTACGGCGATCGATCGATTACTCGCT",
            alignment={
                2: Edit(ins=[Insertion(7, "AA")], subs=[], dels=[Deletion(13, 3)]),
                3: Edit(ins=[Insertion(31, "CCC")], subs=[], dels=[]),
            },
        ),
        2: Block(
            id=1,
            consensus="CATGCTACGCTACGCATTATCGATCGCATCGA",
            alignment={
                5: Edit(ins=[], subs=[], dels=[Deletion(13, 3)]),
                6: Edit(ins=[Insertion(32, "AAA")], subs=[], dels=[]),
            },
        ),
        3: Block(
            id=1,
            consensus="ACGTGTACTAGTACTGC",
            alignment={
                8: Edit(ins=[], subs=[Substitution(12, "C")], dels=[]),
            },
        ),
    }


def test_simplify(graph, expected_graph):
    simplify(graph, focal_paths=["pathA", "pathB"])

    assert graph.paths == expected_graph.paths
    assert graph.blocks == expected_graph.blocks
