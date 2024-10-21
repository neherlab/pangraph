import numpy as np
from . import class_alignments as pga


class Block:
    """Python wrapper for pangraph block object. It has attributes:
    - id (str): block id
    - sequence (str): sequence of the block
    - gaps (list): list of gap occurrences in the alignment, as returned by
        pangraph
    - alignment (pan_alignment): object containing the information provided by
        pangraph that can be used to build alignments. See `pan_alignment`
        class for details.
    """

    def __init__(self, pan_block):
        self.id = pan_block["id"]
        self.sequence = pan_block["sequence"]
        self.alignment = pga.pan_alignment(pan_block)

    def __len__(self):
        """Length of the sequence in base-pairs."""
        return len(self.sequence)

    def __str__(self):
        return f"block {self.id}, consensus len {len(self.sequence)/1000} kbp, {self.depth()} occurrences."

    def depth(self):
        """How many occurrences of the block are present"""
        return len(self.alignment.occs)

    def frequency(self):
        """In how many unique strains is the block present"""
        return len(self.strains())

    def strains(self):
        """Which (unique) strains have the block."""
        return np.unique([w[0] for w in self.alignment.occs])

    def is_duplicated(self):
        """Check whether a block is duplicated in any strain"""
        return self.depth() > self.frequency()

    def generate_alignment(self):
        return self.alignment.to_biopython_aln()
