# Class containing information on a pangraph block alignment.
# It contains utilities to reconstruct the sequences and their alignment

from dataclasses import dataclass
from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import numpy as np


@dataclass
class Substitution:
    pos: int
    alt: str


@dataclass
class Insertion:
    pos: int
    seq: str


@dataclass
class Deletion:
    pos: int
    length: int


class Edits:
    """This class stores edits over the consensus (mutations, insertions, deletions)
    for a block. It can be used to reconstruct the alignment of the block."""

    def __init__(
        self, subs: list[Substitution], ins: list[Insertion], dels: list[Deletion]
    ):
        self.subs = subs
        self.inss = ins
        self.dels = dels

    @staticmethod
    def from_dict(edits: dict) -> "Edits":
        subs = sorted(
            [Substitution(s["pos"], s["alt"]) for s in edits["subs"]],
            key=lambda x: x.pos,
        )
        inss = sorted(
            [Insertion(i["pos"], i["seq"]) for i in edits["inss"]], key=lambda x: x.pos
        )
        dels = sorted(
            [Deletion(d["pos"], d["len"]) for d in edits["dels"]], key=lambda x: x.pos
        )
        return Edits(subs, inss, dels)

    def __str__(self):
        return f"Subs: {self.subs}\nInss: {self.inss}\nDels: {self.dels}"

    def apply_edits(self, cons: str) -> str:
        """Applies the edits to the consensus sequence and returns the edited sequence"""
        sq = list(cons)
        for S in self.subs:
            sq[S.pos] = S.alt
        for D in self.dels:
            for idx in range(D.length):
                sq[D.pos + idx] = ""
        for I in self.inss:
            if I.pos > 0:
                sq[I.pos - 1] += I.seq
            elif I.pos == 0:
                sq[0] = I.seq + sq[0]
        return "".join(sq)

    def aligned_seq(self, cons: str) -> str:
        """Returns the sequence for the alignment, without insertions"""
        sq = list(cons)
        for S in self.subs:
            sq[S.pos] = S.alt
        for D in self.dels:
            for idx in range(D.length):
                sq[D.pos + idx] = "-"
        return "".join(sq)


class Alignment:
    """This class contains information on the sequences contained in a block.
    It saves all the edits between them (mutations, insertions, deletions),
    and can be used to reconstruct the alignment.
    """

    def __init__(self, consensus: str, edits: dict[int, Edits]):
        self.consensus = consensus
        self.edits = edits

    @staticmethod
    def from_dict(block: dict) -> "Alignment":
        consensus = block["consensus"]
        edits = {
            node_id: Edits.from_dict(e) for node_id, e in block["alignments"].items()
        }
        return Alignment(consensus, edits)

    def __len__(self):
        """Returns the number of sequences in the block"""
        return len(self.edits)

    def __str__(self):
        return f"Alignment with {len(self)} sequences and consensus length {len(self.consensus)} bp"

    def __repr__(self):
        return f"Alignment with {len(self)} sequences and consensus length {len(self.consensus)} bp"

    def node_ids(self):
        """Returns the list of node ids"""
        return list(self.edits.keys())

    def generate_sequences(self):
        """returns the unaligned set of reconstructed sequenes for the block"""
        return {
            node_id: e.apply_edits(self.consensus) for node_id, e in self.edits.items()
        }

    def generate_alignment(self) -> dict[int, str]:
        """returns a dictionary node_id -> aligned sequence, excluding insertions."""
        return {
            node_id: e.aligned_seq(self.consensus) for node_id, e in self.edits.items()
        }

    def depth(self):
        """Returns the number of occurrences of the block"""
        return len(self.edits)

    def to_biopython_alignment(self):
        """Returns the alignment in biopython MultipleSeqAlignment format"""
        records = []
        for node_id, seq in self.generate_alignment().items():
            records.append(SeqRecord(Seq(seq), id=str(node_id), description=""))
        return MultipleSeqAlignment(records)

    def to_biopython_records(self):
        """Returns the sequences in biopython SeqRecord format"""
        records = []
        for node_id, seq in self.generate_sequences().items():
            records.append(SeqRecord(Seq(seq), id=str(node_id), description=""))
        return records
