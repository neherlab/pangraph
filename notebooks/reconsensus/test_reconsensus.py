import pytest
from copy import deepcopy

from utils import Edit, Insertion, Deletion, Substitution

from reconsensus import (
    Block,
    Node,
    Path,
    Pangraph,
    find_empty_nodes,
    remove_nodes_from_graph,
    reconsensus,
    reconsensus_mutations,
    reconstruct_deletions,
    reconsensus_insertions,
)


@pytest.fixture
def mut_test_case():
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
def apply_mutations_reconsensus(mut_test_case):
    graph = deepcopy(mut_test_case)
    reconsensus_mutations(graph)
    return graph


class TestReconsensusMut:
    def test_consensus(self, apply_mutations_reconsensus):
        block = apply_mutations_reconsensus
        assert block.consensus == "ATGACTTCGATCTATTCAGAGAA"

    def test_muts_1(self, apply_mutations_reconsensus):
        ed1 = apply_mutations_reconsensus.alignment[1]
        assert ed1.subs == []

    def test_muts_2(self, apply_mutations_reconsensus):
        ed2 = apply_mutations_reconsensus.alignment[2]
        assert ed2.subs == [
            Substitution(pos=10, alt="C"),
            Substitution(pos=17, alt="G"),
        ]

    def test_muts_3(self, apply_mutations_reconsensus):
        ed3 = apply_mutations_reconsensus.alignment[3]
        assert ed3.subs == [
            Substitution(pos=10, alt="C"),
            Substitution(pos=17, alt="G"),
        ]

    def test_muts_4(self, apply_mutations_reconsensus):
        ed4 = apply_mutations_reconsensus.alignment[4]
        assert ed4.subs == [Substitution(pos=1, alt="C")]

    def test_muts_5(self, apply_mutations_reconsensus):
        ed5 = apply_mutations_reconsensus.alignment[5]
        assert ed5.subs == [Substitution(pos=1, alt="G")]

    @pytest.mark.parametrize("idx", range(1, 6))
    def test_no_changes_ins(self, apply_mutations_reconsensus, mut_test_case, idx):
        ed_old = mut_test_case.alignment[idx]
        ed_new = apply_mutations_reconsensus.alignment[idx]
        assert ed_old.ins == ed_new.ins

    @pytest.mark.parametrize("idx", range(1, 6))
    def test_no_changes_dels(self, apply_mutations_reconsensus, mut_test_case, idx):
        ed_old = mut_test_case.alignment[idx]
        ed_new = apply_mutations_reconsensus.alignment[idx]
        assert ed_old.dels == ed_new.dels
