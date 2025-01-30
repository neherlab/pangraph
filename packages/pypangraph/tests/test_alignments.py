import pytest
import pypangraph as pp
import numpy as np
from pypangraph.class_alignments import (
    Alignment,
    Edits,
    Insertion,
    Deletion,
    Substitution,
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


def test_reconstruct_sequences(alignment):
    seqs = alignment.generate_sequences()
    assert seqs[1] == "AGTCTACCTACTGGCAC"
    assert seqs[2] == "ACTCTACCCGCTACTAGCACAAA"
    assert seqs[3] == "ACTGGCTACCCGACTGGCAC"


def test_reconstruct_alignment(alignment):
    aln = alignment.generate_alignment()
    assert aln[1] == "AGTCTACC---TACTGGCAC"
    assert aln[2] == "ACTCTACCCGCTACTAGCAC"
    assert aln[3] == "ACTCTACCCG--ACTGGCAC"


@pytest.fixture
def graph():
    fname = "tests/data/plasmids.json"
    pan = pp.Pangraph.from_json(fname)
    return pan


def test_get_sequences(graph):
    bid = graph.blocks.keys()[0]
    block = graph.blocks[bid]
    seqs = block.to_sequences()
    assert len(seqs) == 15


def test_get_alignment(graph):
    bid = graph.blocks.keys()[0]
    block = graph.blocks[bid]
    aln = block.to_alignment()
    # all sequences have the same length
    assert len(set(len(seq) for seq in aln.values())) == 1
    # all sequences have the same length as the consensus
    assert len(block.consensus()) == len(list(aln.values())[0])


def test_core_alignment(graph):
    aln = graph.core_genome_alignment()
    A = np.array(aln)
    assert A.shape[0] == 15
    assert A.shape[1] == 64989
