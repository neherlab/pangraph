import pytest
from copy import deepcopy

from utils import (
    Edit,
    Insertion,
    Deletion,
    Substitution,
    map_variations,
)

from reconsensus import (
    Block,
    Node,
    Path,
    Pangraph,
    apply_indels,
    majority_deletions,
    majority_insertions,
    reconsensus_mutations,
    reconsensus,
)


@pytest.fixture
def block_1():
    consensus = "AGGACTTCGATCTATTCGGAGAA"
    #            0         1         2
    #            01234567890123456789012
    # node 1)    .T...--..........A.....
    # node 2)    .T...--...C......|.....
    # node 3)    .T...--...C............
    # node 4)    .C.......---.....A.....
    # node 5)    .....|...........A.....
    #  L = 23, N = 5

    aln = {
        1: Edit(
            subs=[Substitution(pos=1, alt="T"), Substitution(pos=17, alt="A")],
            dels=[Deletion(pos=5, length=2)],
            ins=[],
        ),
        2: Edit(
            subs=[Substitution(pos=1, alt="T"), Substitution(pos=10, alt="C")],
            dels=[Deletion(pos=5, length=2)],
            ins=[Insertion(pos=17, ins="CAAT")],
        ),
        3: Edit(
            subs=[Substitution(pos=1, alt="T"), Substitution(pos=10, alt="C")],
            dels=[Deletion(pos=5, length=2)],
            ins=[],
        ),
        4: Edit(
            subs=[Substitution(pos=1, alt="C"), Substitution(pos=17, alt="A")],
            dels=[Deletion(pos=9, length=3)],
            ins=[],
        ),
        5: Edit(
            subs=[Substitution(pos=17, alt="A")],
            dels=[Deletion(pos=5, length=2)],
            ins=[],
        ),
    }

    block = Block(0, consensus, aln)
    return block


@pytest.fixture
def block_1_mut_reconsensus():
    consensus = "ATGACTTCGATCTATTCAGAGAA"
    #            0         1         2
    #            01234567890123456789012
    # node 1)    .....--................
    # node 2)    .....--...C......G|....
    # node 3)    .....--...C......G.....
    # node 4)    .C.......---...........
    # node 5)    .G...|.................
    #  L = 23, N = 5

    aln = {
        1: Edit(
            subs=[],
            dels=[Deletion(pos=5, length=2)],
            ins=[],
        ),
        2: Edit(
            subs=[Substitution(pos=10, alt="C"), Substitution(pos=17, alt="G")],
            dels=[Deletion(pos=5, length=2)],
            ins=[Insertion(pos=17, ins="CAAT")],
        ),
        3: Edit(
            subs=[Substitution(pos=10, alt="C"), Substitution(pos=17, alt="G")],
            dels=[Deletion(pos=5, length=2)],
            ins=[],
        ),
        4: Edit(
            subs=[Substitution(pos=1, alt="C")],
            dels=[Deletion(pos=9, length=3)],
            ins=[],
        ),
        5: Edit(
            subs=[Substitution(pos=1, alt="G")],
            dels=[Deletion(pos=5, length=2)],
            ins=[],
        ),
    }

    block = Block(0, consensus, aln)
    return block


@pytest.fixture
def block_2():
    consensus = "AGGACTTCGATCTATTCGGAGAA"
    #            0         1         2
    #            01234567890123456789012
    # node 1)   |.T.|.--......|...A..-..
    # node 2)   |.T...--...C..|......--.|
    # node 3)    .T..----..C............|
    # node 4)    .C.|.....---.....A.....|
    # node 5)   |...|.........|...A.--..
    #  L = 23, N = 5

    aln = {
        1: Edit(
            subs=[Substitution(1, "T"), Substitution(17, "A")],
            dels=[Deletion(5, 2), Deletion(20, 1)],
            ins=[Insertion(0, "G"), Insertion(3, "AA"), Insertion(13, "AA")],
        ),
        2: Edit(
            subs=[Substitution(1, "T"), Substitution(10, "C")],
            dels=[Deletion(5, 2), Deletion(20, 2)],
            ins=[Insertion(0, "G"), Insertion(13, "AA"), Insertion(23, "TT")],
        ),
        3: Edit(
            subs=[Substitution(1, "T"), Substitution(10, "C")],
            dels=[Deletion(4, 4)],
            ins=[Insertion(23, "TT")],
        ),
        4: Edit(
            subs=[Substitution(1, "C"), Substitution(17, "A")],
            dels=[Deletion(9, 3)],
            ins=[Insertion(3, "C"), Insertion(23, "TT")],
        ),
        5: Edit(
            subs=[Substitution(17, "A")],
            dels=[Deletion(19, 2)],
            ins=[Insertion(3, "C"), Insertion(0, "G"), Insertion(13, "AA")],
        ),
    }

    block = Block(0, consensus, aln)
    return block


def test_reconsensus_mutations(block_1, block_1_mut_reconsensus):
    reconsensus_mutations(block_1)
    assert block_1 == block_1_mut_reconsensus


def test_majority_deletions(block_2):
    dels = majority_deletions(block_2)
    assert dels == [5, 6, 20]


def test_majority_insertions(block_2):
    ins = majority_insertions(block_2)
    assert ins == [(0, "G"), (13, "AA"), (23, "TT")]


def test_apply_indels(block_2):
    dels = majority_deletions(block_2)
    ins = majority_insertions(block_2)
    cons = block_2.consensus
    cons = apply_indels(cons, dels, ins)
    assert cons == "GAGGACCGATCTAAATTCGGAAATT"


def test_map_variations():
    """
    NB: No need to reproduce this test, just for me to test that my
    re-implementation works
    """
    consensus = "AGGACTTCGATCTATTCGG"
    qry = "AGGACCGATCTAAATACGG"
    edits = map_variations(consensus, qry)
    assert edits == Edit(
        ins=[Insertion(13, "AA")],
        dels=[Deletion(5, 2)],
        subs=[Substitution(15, "A")],
    )


@pytest.fixture
def block_3():
    consensus = "GCCTCTTCCCGACCACGCGTTACAACATGGGACAGGCCTGCGCTTGAGGC"
    #            0         1         2         3         4
    #            01234567890123456789012345678901234567890123456789
    # node 1)   ......A.............----...........................
    # node 2)   ......A..............---............|..............|
    # node 3)   ...............G............G......................
    # node 4)   ......A..............---...........................|
    # node 5)   ...................................................|
    #  L = 50, N = 5
    aln = {
        1: Edit(
            subs=[Substitution(5, "A")],
            dels=[Deletion(19, 4)],
            ins=[],
        ),
        2: Edit(
            subs=[Substitution(5, "A")],
            dels=[Deletion(20, 3)],
            ins=[Insertion(35, "AA"), Insertion(50, "TT")],
        ),
        3: Edit(
            subs=[Substitution(14, "G"), Substitution(27, "G")],
            dels=[],
            ins=[],
        ),
        4: Edit(
            subs=[Substitution(5, "A")],
            dels=[Deletion(20, 3)],
            ins=[Insertion(50, "TT")],
        ),
        5: Edit(
            subs=[],
            dels=[],
            ins=[Insertion(50, "TT")],
        ),
    }
    return Block(3, consensus, alignment=aln)


@pytest.fixture
def block_3_reconsensus():
    consensus = "GCCTCATCCCGACCACGCGTAACATGGGACAGGCCTGCGCTTGAGGCTT"
    #            0         1         2         3         4
    #            0123456789012345678901234567890123456789012345678
    # node 1)   ....................-...........................--
    # node 2)   .................................|................
    # node 3)   ......T........G.....|..........................--
    # node 4)   ..................................................
    # node 5)   ......T..............|............................
    #  L = 50, N = 5
    aln = {
        1: Edit(
            subs=[],
            dels=[Deletion(19, 1), Deletion(47, 2)],
            ins=[],
        ),
        2: Edit(
            subs=[],
            dels=[],
            ins=[Insertion(32, "AA")],
        ),
        3: Edit(
            subs=[Substitution(5, "T"), Substitution(14, "G"), Substitution(24, "G")],
            dels=[Deletion(47, 2)],
            ins=[Insertion(20, "TAC")],
        ),
        4: Edit(
            subs=[],
            dels=[],
            ins=[],
        ),
        5: Edit(
            subs=[Substitution(5, "T")],
            dels=[],
            ins=[Insertion(20, "TAC")],
        ),
    }
    return Block(3, consensus, alignment=aln)


def test_reconsensus(block_3, block_3_reconsensus):
    reconsensus(block_3)
    assert block_3.consensus == block_3_reconsensus.consensus
    assert block_3.alignment == block_3_reconsensus.alignment
