# Utilities to deal with biopython objects, such as sequences, alignments and
# phylogenetic trees

import numpy as np
import matplotlib as mpl

from Bio import SeqIO
from Bio import Seq
from Bio import AlignIO

# ------------------- Alignments ----------------------------


def seqmat_to_align(seqMat, seqIds):
    """
    given a numpy array of character of shape NxL, containing N sequences, it
    turns it into an alignment object, whose sequence ids are specified by the
    second argument
    """
    seqRecords = []
    for n, x in enumerate(seqMat):
        seq = Seq.Seq("".join(x))
        rec = SeqIO.SeqRecord(seq, id=seqIds[n])
        seqRecords.append(rec)
    alignment = AlignIO.MultipleSeqAlignment(seqRecords)
    return alignment


# ------------------- Phylogenetic Trees ----------------------------


def initialize_property(tree, prop_name, ivalue):
    """given a tree from biopython's Phylo module, it initialzies for each of
    the clades a property with the specified name, and initializes it to the
    specified value
    """

    def rec_initialize(node):
        setattr(node, prop_name, ivalue)
        for child in node:
            rec_initialize(child)

    rec_initialize(tree.root)


def max_property(tree, prop_name, exclude_terminal=False):
    """Finds the maximum value of a property on the tree. Optionally,
    terminal nodes can be excluded
    """

    def rec_max(node):
        values = [getattr(node, prop_name)]
        for child in node:
            if (not exclude_terminal) or (not child.is_terminal()):
                values.append(rec_max(child))
        return np.max(values)

    return rec_max(tree.root)


def tot_property(tree, prop_name, exclude_terminal=False):
    """Finds the sum of values of a property on the tree. Optionally,
    terminal nodes can be excluded
    """

    def rec_tot(node):
        values = [getattr(node, prop_name)]
        for child in node:
            if (not exclude_terminal) or (not child.is_terminal()):
                values.append(rec_tot(child))
        return np.sum(values)

    return rec_tot(tree.root)


def color_property(tree, prop_name, colorf, ignore_terminal=False):
    """sets the `color` attribute of clades based on the specified clade property
    and the specified color function colorf(property) -> color.
    If ignore_terminal is set to True then terminal nodes are colored in gray.
    """

    def rec_color(node):
        # if terminal leaves should be ignored
        if ignore_terminal:
            if node.is_terminal():
                node.color = mpl.colors.to_hex("lightgray")
                return

        prop_v = getattr(node, prop_name)
        color = colorf(prop_v)
        node.color = mpl.colors.to_hex(color)
        for child in node:
            color_property(child, prop_name, colorf, ignore_terminal)

    rec_color(tree.root)


def is_monophyletic(tree, leaves_ids):
    """given a list of leaves ids, checks whether there is a single clade whose
    terminal nodes are exactly the specified leaves ids (i.e. if they are
    monophyletic).
    Returns False if such a clade was not found, and otherwise it returns the
    clade.
    """

    # set of all leaves
    leaves = tree.get_terminals()

    # subsets of leaves considered
    lvs = [l for l in leaves if l.name in leaves_ids]
    cl = tree.is_monophyletic(lvs)
    return cl
