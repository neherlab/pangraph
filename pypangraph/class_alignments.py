# Class containing information on a pangraph block alignment.
# It contains utilities to reconstruct the alignment, the sequences
# And the set of SNPs in the columns of the alignment without gaps

import numpy as np
from Bio import SeqRecord, Seq, AlignIO


class pan_alignment:
    """This class contains information on the sequences contained in a block.
    It saves all the variation between them (mutations, insertions, deletions),
    and can be used to reconstruct the alignment.

    It has attributes:
    - consensus: the block consensus sequence
    - occs: a list (strain, n. occurrence, strand) of block occurrences
    - gaps: a list [[pos, len], ...] of gaps to be added to the consensus
    - muts: a dictionary {occurrence -> mutations}. Each mutation is saved as
        [position, nucleotide]
    - ins: a dictionary {occurrence -> insertions}. Each insertion is saved as
        [[sequence position, gap position], sequence]. Insertions are to be inserted
        in a gap. The first position indicates the gap position on the consensus,
        and the second position indicates the beginning od the insertion w.r.t. the
        gap (0 means at the beginning of the gap).
    - dels: a dictionary {occurrence -> deletions}. Each deletion is saved as
        [position, deletion length]
    - pos: a dictionary {occurrence -> [beg, end]}, where `beg` and `end` are julia
        indices indicating where the block occurrence lays on the full genome
        sequence. in python [i:j] -> [i-1:j].
        NB: for blocks that wrap around the end of the genome, it can be j<i!
    """

    def __init__(self, pan_block: dict):
        self.consensus = pan_block["sequence"]
        self.gaps = pan_block["gaps"]
        occs, muts, ins, dels, pos = parse_alignment(pan_block)
        self.occs = occs
        self.muts = muts
        self.ins = ins
        self.dels = dels
        self.pos = pos  # [i:j] as per julia numbering, to python is [i-1:j]

    def block_occurrence_length(self, occ: tuple):
        """Returns the length of a particular block occurrence, inferred from alignment
        information."""
        L = len(self.consensus)
        ins, dels = self.ins[occ], self.dels[occ]
        for i in ins:
            L += len(i[1])  # add insertions
        for d in dels:
            L -= d[1]  # remove deletions
        return L

    def generate_alignments(self, which=None):
        """Returns the aligned set of sequences corresponding to the same block,
        together with the corresponding list of occurrences (strain, occurrence_n, strand).
        Optionally a subset of occurrences can be selected. These must be elements of `self.occs`
        """
        if which is None:
            which = self.occs

        seqs = []
        for wh in which:
            seq = reconstruct_alignment(
                self.consensus,
                gaps=self.gaps,
                muts=self.muts[wh],
                ins=self.ins[wh],
                dels=self.dels[wh],
            )
            seqs.append(seq)
        return seqs, which

    def generate_sequences(self, which=None):
        """Returns the non-aligned set of sequences corresponding to the same block,
        together with the corresponding list of occurrences (strain, occurrence_n, strand).
        Optionally a subset of occurrences can be selected. These must be elements of `self.occs`
        """
        alns, which = self.generate_alignments(which)
        seqs = [aln.replace("-", "") for aln in alns]  # remove gaps "-"
        return seqs, which

    def extract_nongap_SNPs(self, which=None):
        """Given a set of genomes, it returns the positions in which they have mutations w.r.t the block
        consensus.

        Returns:
            - positions (np.array): list of positions in which mutations occurr (nb: relative to
                block consensus)
            - SNPs (np matrix): matrix of nucleotides. The first index represent occurrences and
                the second the position
            - which (list): set of occurrences corresponding to the matrix rows.

        NB: some alignment columns might be the same if a subset of strains is selected. In this case
            a warning message is printed.
        NB: the order of the occurrences can be controlled through the optional `which argument`.
        """

        if which is None:
            which = self.occs

        positions, SNPs = extract_relevant_SNPs(
            consensus=self.consensus, dels=self.dels, muts=self.muts, which=which
        )

        return positions, SNPs, which

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


def parse_alignment(pan_block: dict):
    """Given the pan_block dictionary, parses the results and returns the list of
    - occurrences (block_name, occurrence_number, strans)
    - mutations
    - insertions
    - deletions
    - position
    """

    # define containers
    pan_labels = ["mutate", "insert", "delete", "positions"]
    muts, ins, dels, pos = {}, {}, {}, {}
    containers = [muts, ins, dels, pos]

    # capture occurrences
    occs = []
    for occ, _ in pan_block["mutate"]:
        occ_id = blockid_to_tuple(occ)
        occs.append(occ_id)

    for pan_label, container in zip(pan_labels, containers):
        for occ, item in pan_block[pan_label]:
            occ_id = blockid_to_tuple(occ)
            assert occ_id in occs, "mismatch in occurrences"
            container[occ_id] = item

    return occs, muts, ins, dels, pos


def blockid_to_tuple(bl):
    """Given a block id, returns a tuple that can be used as key in a dictionary"""
    return (bl["name"], bl["number"], bl["strand"])


def reconstruct_alignment(
    consensus: str, gaps: dict, muts: list, ins: list, dels: list
) -> str:
    """Given a consensus sequence, it applies all the mutations/insertions/deletions to
    reconstruct the alignment, and returns the sequence + gaps"""
    seq = list(consensus)

    if len(seq) == 0:
        seq = [""]
        print("WARNING: block of length zero present")

    # apply mutations
    for pos, nt in muts:
        seq[pos - 1] = nt

    # apply deletions (insert gaps)
    for pos, L in dels:
        for l in range(L):
            seq[pos - 1 + l] = "-"

    # create dictionary of gaps, to later be inserted
    gap_dict = {}
    for pos, L in gaps.items():
        gap_dict[int(pos)] = ["-"] * L

    # fill these gaps with insertions
    for ins_descr, nts in ins:
        gap_id, gap_pos = ins_descr
        gap = gap_dict[gap_id]  # capture gap
        for i, nt in enumerate(nts):
            gap[gap_pos + i] = nt

    # add the gaps to the sequence
    for i, gap in gap_dict.items():
        gap_nt = "".join(gap)
        if i > 0:  # add after position i (julia) or after i-1 (python)
            seq[i - 1] = seq[i - 1] + gap_nt
        else:  # add at the beginning
            seq[0] = gap_nt + seq[0]

    # concatenate and return the sequence
    seq = "".join(seq)
    return seq


def extract_relevant_SNPs(consensus: str, dels: list, muts: list, which: list):
    """Function to extract a set of relevant SNPs from the set of mutations and deletions.
    Optionally a particular subset of strains can be considered. It only takes columns of
    the alignment without deletions. It returns:
    - positions (np.array): list of positions in which mutations occurr (nb: relative to
        block consensus, and in python indexing)
    - SNPs (np matrix): matrix of nucleotides. The first index represent occurrences and
        the second the position
    - which (list): set of occurrences corresponding to the matrix rows.

    Nb: columns might still be equal to consensus, but only if a subset of strains (all
        having the same mutations) is requested.
    Nb: occurrences are in the same order as passed in the `which argument`
    Nb: the position is relative to the BLOCK consensus. this is not the same as the
        alignment consensus.
    """

    # find all positions with a deletions. These cannot be considered for
    # valid SNPs since they would contain a gap in one column
    forbidden_positions = set()  # julia indexing
    for wh in which:
        ds = dels[wh]
        for pos, L in ds:
            forbidden_positions |= set([pos + l for l in range(L)])

    # find relevant columns
    snp_pos = []
    snp_columns = {}
    W = len(which)
    for i, wh in enumerate(which):
        mt = muts[wh]
        for pos, nt in mt:
            # exclude positions that have gaps
            if pos in forbidden_positions:
                continue

            if pos in snp_pos:
                # if position was already present, add the mutation in the specific column
                snp_columns[pos][i] = nt
            else:
                # else add the new position and add a new column
                snp_pos.append(pos)
                snp_columns[pos] = [
                    consensus[pos - 1]
                ] * W  # NB: pos is in julia indexing
                snp_columns[pos][i] = nt

    # format and return relevant results
    positions = np.sort(snp_pos)
    SNPs = np.zeros((W, positions.size), dtype=str)
    for i, pos in enumerate(positions):
        SNPs[:, i] = snp_columns[pos]

    # print a warning if there is at least one column where all are consensus
    is_there_a_consensus = np.any(np.all(SNPs == SNPs[0, :], axis=0))  # TODO: test?
    if is_there_a_consensus:
        print("Warning: at least one column has no mutations")
        print("Maybe it corresponds to a subset?")

    positions -= 1  # julia to python indexing
    return positions, SNPs
