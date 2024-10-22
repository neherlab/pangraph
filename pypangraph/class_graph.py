# Wrapper class to load a Pangraph object from a .json file.

import numpy as np
import json
import pandas as pd
import jsonschema
from Bio import SeqRecord, Seq, AlignIO

from collections import Counter, defaultdict

from .class_path import PathCollection
from .class_block import BlockCollection
from .class_node import Nodes
from .pangraph_schema import schema


class Pangraph:
    """Wrapper class to load and interact with the output of the Pangraph pipeline.
    The class has three main attributes:
    - `paths` : each strain has a path representaiton. A path consists in a series of blocks.
    - `blocks` : the alignment of each homologous sequence is represented as a block.
    - `nodes` : each node represents a particular occurrence of a block on a path.
    """

    def __init__(self, pan_json):
        """Python calss to load the output of the Pangraph pipeline.

        Args:
            pan_json (dict): content of the .json file produced by pangraph.
        """
        self.paths = PathCollection(pan_json["paths"])
        self.blocks = BlockCollection(pan_json["blocks"])
        self.nodes = Nodes(pan_json["nodes"])

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
        return list(self.paths.keys())

    def to_blockcount_df(self):
        """Returns a dataframe whose rows are strain names, and columns are block
        names. Values indicate the number of times a block is present. This can
        also be used to build a presence / absence matrix."""
        df = self.nodes.block_path_counts()
        path_name_dict = {path.id: path.name for name, path in self.paths.items()}
        df.rename(columns=path_name_dict, inplace=True)
        return df

    def to_blockstats_df(self):
        """Returns a dataframe containing statistics about blocks distribution.
        The index of the dataframe are block ids, and the columns are:
        - count: n. times the block occurs
        - n. strains: number of strains in which the block is observed
        - duplicated: whether the block is duplicated in at least one strain
        - len: average block length from pangraph.
        - core: whether a gene occurrs exactly once per strain
        """
        df = self.nodes.to_blockstats_df()
        df["len"] = [len(self.blocks[bid]) for bid in df.index]
        return df

    def core_genome_alignment(self, guide_strain=None):
        # TODO
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
