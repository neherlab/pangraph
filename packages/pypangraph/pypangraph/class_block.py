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

    def __init__(self, block_id: int, alignment: pga.Alignment):
        self.id = block_id
        self.alignment = alignment

    @staticmethod
    def from_dict(block: dict) -> "Block":
        block_id = block["id"]
        alignment = pga.Alignment.from_dict(block)
        return Block(block_id, alignment)

    def __len__(self):
        """Length of the sequence in base-pairs."""
        return len(self.alignment.consensus)

    def __str__(self):
        return f"block {self.id}, consensus len = {len(self.alignment.consensus)} bp, n. nodes = {self.depth()}"

    def __repr__(self):
        return self.__str__()

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

    def to_biopython_alignment(self):
        """Returns the block alignment as a Biopython alignment object"""
        return self.alignment.to_biopython_alignment()

    def to_biopython_records(self):
        """Returns the block sequence as a list of Biopython SeqRecord objects"""
        return self.alignment.to_biopython_records()


class BlockCollection(IndexedCollection):
    """Collection of all blocks. Inherits from IndexedCollection to allow for
    smart indexing of blocks.
    """

    def __init__(self, pan_blocks):
        ids = [block["id"] for block in pan_blocks.values()]
        items = [Block.from_dict(block) for block in pan_blocks.values()]
        IndexedCollection.__init__(self, ids, items)
