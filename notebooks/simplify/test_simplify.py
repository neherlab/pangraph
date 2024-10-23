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
    #   n1:    ..........G.....................  l = 32
    #   n2:    .............xxx................  l = 29
    #   n3:    ................................  l = 32
    return Block(
        id=1,
        consensus="CATGCTACGCTACGCATTATCGATCGCATCGA",
        alignment={
            1: Edit(ins=[], subs=[Substitution(10, "G")], dels=[]),
            2: Edit(ins=[], subs=[], dels=[Deletion(13, 3)]),
            3: Edit(ins=[Insertion(32, "CCC")], subs=[], dels=[]),
        },
    )
