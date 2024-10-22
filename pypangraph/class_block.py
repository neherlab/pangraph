import numpy as np
from . import class_alignments as pga


class Block:
    """Pangraph block object. It has attributes:
    - id (str): block id
    - sequence (str): consensus sequence of the block
    - alignment (pan_alignment): object containing the information provided by
        pangraph that can be used to build alignments. See `pan_alignment`
        class for details.
    """

    def __init__(self, block):
        self.id = block["id"]
        self.sequence = block["consensus"]
        self.alignment = pga.Alignment(block)

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
