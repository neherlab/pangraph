import numpy as np
from . import class_alignments as pga
from .indexed_collection import IndexedCollection


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
        self.alignment = pga.Alignment(block)

    def __len__(self):
        """Length of the sequence in base-pairs."""
        return len(self.alignment.consensus)

    def __str__(self):
        return f"block {self.id}, consensus len {len(self.sequence)/1000} kbp, {self.depth()} occurrences."

    def depth(self):
        """How many occurrences of the block are present"""
        return self.alignment.depth()

    def consensus(self) -> str:
        """Returns the consensus sequence of the block"""
        return self.alignment.consensus

    def to_sequences(self) -> dict[int, str]:
        """Returns a dictionary node_id -> sequence for the block"""
        return self.alignment.generate_sequences()

    def to_alignment(self) -> dict[int, str]:
        """Returns a dictionary node_id -> aligned sequence for the block.
        The aligned sequence does not include insertions."""
        return self.alignment.generate_alignment()


class BlockCollection(IndexedCollection):
    """Collection of all blocks. Inherits from IndexedCollection to allow for
    smart indexing of blocks.
    """

    def __init__(self, pan_blocks):
        ids = [block["id"] for block in pan_blocks.values()]
        items = [Block(block) for block in pan_blocks.values()]
        IndexedCollection.__init__(self, ids, items)
