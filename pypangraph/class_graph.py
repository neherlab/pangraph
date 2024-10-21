# Wrapper to import in python the results of the Pangraph pipeline.

import numpy as np
import json
import pandas as pd
import jsonschema
from Bio import SeqRecord, Seq, AlignIO

from collections import Counter, defaultdict

from .indexed_collection import PathCollection, BlockCollection, NodeCollection
from .pangraph_schema import schema


class Pangraph:
    """Wrapper class to load and interact with the output of the Pangraph pipeline.
    The class has two main attributes:
    - `paths` : each strain has a path representaiton. A path consists in a series of blocks.
    - `blocks` : the alignment of each homologous sequence is represented as a block.
    """

    def __init__(self, pan_json):
        """Python calss to load the output of the Pangraph pipeline.

        Args:
            pan_json (dict): content of the .json file produced by pangraph.
        """
        self.paths = PathCollection(pan_json["paths"])
        self.blocks = BlockCollection(pan_json["blocks"])
        self.nodes = NodeCollection(pan_json["nodes"])

    @staticmethod
    def load_json(filename):
        """Creates a Pangraph object by loading it from the .json file.

        Args:
            load_json (str): .json file to be loaded.

        Returns:
            Pangraph: the Pangraph object containing the results of the pipeline.
        """

        isjson = str(filename).endswith(".json")
        if not isjson:
            raise Exception(f"the input file {filename} should be in .json format")

        with open(filename, "r") as f:
            pan_json = json.load(f)

        try:
            graph = {"pangraph": pan_json}
            jsonschema.validate(instance=graph, schema=schema)
        except jsonschema.exceptions.ValidationError as ex:
            print(ex)

        pan = Pangraph(pan_json)
        return pan

    def strains(self):
        """Return lists of strain names"""
        return self.paths.ids_copy()

    def block_ids(self):
        """Returns the list of block ids"""
        return self.blocks.ids_copy()

    def to_paths_dict(self):
        """Generates a compressed representation of paths as simply lists of
        block ids. Returns a dictionary strain name -> list of block ids.
        """
        return self.paths.to_block_dict()

    def to_blockcount_df(self):
        """Returns a dataframe whose rows are strain names, and columns are block
        names. Values indicate the number of times a block is present. This can
        also be used to build a presence / absence matrix."""
        block_counters = {}
        for path in self.paths:
            block_counters[path.name] = Counter(path.block_ids)
        return pd.DataFrame(block_counters).fillna(0).T

    def to_blockstats_df(self):
        """Returns a dataframe containing statistics about blocks distribution.
        The index of the dataframe are block ids, and the columns are:
        - count: n. times the block occurs
        - n. strains: number of strains in which the block is observed
        - duplicated: whether the block is duplicated in at least one strain
        - len: average block length from pangraph.
        - core: whether a gene occurrs exactly once per strain
        """
        block_counter = Counter()
        str_counter = Counter()

        for path in self.paths:
            block_counter.update(path.block_ids)
            str_counter.update(np.unique(path.block_ids))

        lengths = {block.id: len(block.sequence) for block in self.blocks}

        df = pd.DataFrame(
            {"count": block_counter, "n. strains": str_counter, "len": lengths}
        )
        df["duplicated"] = df["count"] > df["n. strains"]
        df["core"] = (df["n. strains"] == len(self.paths)) & (~df["duplicated"])
        return df

    def core_genome_alignment(self, guide_strain=None):
        """Returns the core genome aligment, in a biopython alignment object.
        The order of the blocks is determined by the guide strain, if provided.
        """
        strains = self.strains()

        # get core blocks
        bdf = self.to_blockstats_df()
        core_blocks = bdf[bdf["core"]].index

        # get order of core blocks in guide strain
        if guide_strain is None:
            guide_strain = self.strains()[0]
        assert (
            guide_strain in strains
        ), f"Guide strain {guide_strain} not found in the dataset"
        guide_path = self.paths[guide_strain]
        core_blocks = [
            (bid, strand)
            for bid, strand in zip(guide_path.block_ids, guide_path.block_strands)
            if bid in core_blocks
        ]

        # get alignment
        alignment = defaultdict(str)
        for bid, guide_strand in core_blocks:
            block = self.blocks[bid]
            seqs, which = block.alignment.generate_alignments()
            assert set([w[0] for w in which]) == set(
                strains
            ), f"error: strain missing in block {block}"
            for seq, occ in zip(seqs, which):
                strain, occ_n, strand = occ
                assert (
                    occ_n == 1
                ), f"error in {bid=} | {occ=}, only one occurrence is expected"
                if not guide_strand:
                    seq = str(Seq.Seq(seq).reverse_complement())
                alignment[strain] += seq

        # convert to biopython alignment
        records = []
        for strain, seq in alignment.items():
            record = SeqRecord.SeqRecord(Seq.Seq(seq), id=strain, description="")
            records.append(record)
        alignment = AlignIO.MultipleSeqAlignment(records)

        return alignment
