import pytest
import pypangraph as pp
from pypangraph.class_block import Block
from pypangraph.class_alignments import (
    Alignment,
    Substitution,
    Insertion,
    Deletion,
    Edits,
)


@pytest.fixture
def alignment():
    #        0           1
    #        012  34567890123456789
    # cons   ACT  CTACCCGCTACTGGCAC
    # n1      G        ---
    # n2                      A    |AAA
    # n3        GG       --
    consensus = "ACTCTACCCGCTACTGGCAC"
    n1 = Edits(subs=[Substitution(1, "G")], ins=[], dels=[Deletion(8, 3)])
    n2 = Edits(subs=[Substitution(15, "A")], ins=[Insertion(20, "AAA")], dels=[])
    n3 = Edits(subs=[], ins=[Insertion(3, "GG")], dels=[Deletion(10, 2)])
    return Alignment(consensus, {1: n1, 2: n2, 3: n3})


@pytest.fixture
def block(alignment):
    block_id = 42
    block = Block(block_id, alignment)
    return block


def test_block_stats(block):
    assert len(block) == 20
    assert block.depth() == 3
    assert block.consensus() == "ACTCTACCCGCTACTGGCAC"


def test_block_sequences(block):
    seqs = block.to_sequences()
    assert seqs[1] == "AGTCTACCTACTGGCAC"
    assert seqs[2] == "ACTCTACCCGCTACTAGCACAAA"
    assert seqs[3] == "ACTGGCTACCCGACTGGCAC"


def test_block_alignment(block):
    aln = block.to_alignment()
    assert aln[1] == "AGTCTACC---TACTGGCAC"
    assert aln[2] == "ACTCTACCCGCTACTAGCAC"
    assert aln[3] == "ACTCTACCCG--ACTGGCAC"


def test_block_biopython_alignment(block):
    aln = block.to_biopython_alignment()
    assert len(aln) == 3
    assert aln.get_alignment_length() == 20
    assert aln[0].seq == "AGTCTACC---TACTGGCAC"
    assert aln[0].id == "1"
    assert aln[1].seq == "ACTCTACCCGCTACTAGCAC"
    assert aln[1].id == "2"
    assert aln[2].seq == "ACTCTACCCG--ACTGGCAC"
    assert aln[2].id == "3"


def test_block_biopython_records(block):
    records = block.to_biopython_records()
    assert len(records) == 3
    assert records[0].id == "1"
    assert records[0].seq == "AGTCTACCTACTGGCAC"
    assert records[1].id == "2"
    assert records[1].seq == "ACTCTACCCGCTACTAGCACAAA"
    assert records[2].id == "3"
    assert records[2].seq == "ACTGGCTACCCGACTGGCAC"
