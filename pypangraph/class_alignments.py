# Class containing information on a pangraph block alignment.
# It contains utilities to reconstruct the alignment, the sequences
# And the set of SNPs in the columns of the alignment without gaps

import numpy as np
from Bio import SeqRecord, Seq, AlignIO
from dataclasses import dataclass


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

    def __init__(self, edits: dict) -> None:
        self.subs = sorted(
            [Substitution(s["pos"], s["alt"]) for s in edits["subs"]],
            key=lambda x: x.pos,
        )
        self.inss = sorted(
            [Insertion(i["pos"], i["seq"]) for i in edits["inss"]], key=lambda x: x.pos
        )
        self.dels = sorted(
            [Deletion(d["pos"], d["len"]) for d in edits["dels"]], key=lambda x: x.pos
        )

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
                sq[I.pos - 1] += I.ins
            elif I.pos == 0:
                sq[0] = I.ins + sq[0]
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

    def __init__(self, block: dict):
        self.consensus = block["consensus"]
        self.edits = {node_id: Edits(e) for node_id, e in block["alignments"].items()}

    def node_ids(self):
        """Returns the list of node ids"""
        return list(self.edits.keys())

    def generate_sequences(self, which=None):
        """Returns the non-aligned set of sequences corresponding to the same block,
        together with the corresponding list of occurrences (strain, occurrence_n, strand).
        Optionally a subset of occurrences can be selected. These must be elements of `self.occs`
        """
        alns, which = self.generate_alignments(which)
        seqs = [aln.replace("-", "") for aln in alns]  # remove gaps "-"
        return seqs, which

    def to_biopython_aln(self):
        """Returns a biopython alignment object for the block"""
        seqs, which = self.generate_alignments()
        records = []
        for seq, occ in zip(seqs, which):
            block_id, n, strand = occ
            record = SeqRecord.SeqRecord(
                Seq.Seq(seq), id=f"{block_id=}|{n=}|{strand=}", description=""
            )
            records.append(record)
        return AlignIO.MultipleSeqAlignment(records)
