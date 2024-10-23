import pandas as pd
import numpy as np


def parse_strandedness(strand: str) -> bool:
    """Parse strandedness from a string to a boolean:
    + -> True,
    - -> False
    """
    match strand:
        case "+":
            return True
        case "-":
            return False
        case _:
            raise ValueError(f"Strand {strand} not recognized.")


class Nodes:
    """Dataset with pangraph nodes information:
    - id (str): node id. Index of the dataframe
    - block_id (str): block id
    - path_id (str): path id
    - strand (bool): strandedness of the node
    - position (tuple): position of the node in the genome (start, end)
    """

    def __init__(self, nodes_dict):
        self.df = pd.DataFrame.from_dict(
            {
                node_id: {
                    "block_id": node["block_id"],
                    "path_id": node["path_id"],
                    "strand": parse_strandedness(node["strand"]),
                    "start": node["position"][0],
                    "end": node["position"][1],
                }
                for node_id, node in nodes_dict.items()
            }
        ).T

    def block_path_counts(self) -> pd.DataFrame:
        """Returns a dataframe with columsns=paths, index=blocks, values=counts"""
        return self.df.pivot_table(
            index="block_id", columns="path_id", aggfunc="size", fill_value=0
        )

    def block_counts(self) -> pd.Series:
        """Returns a pandas series with columns=blocks, values=#nodes"""
        return self.df["block_id"].value_counts()

    def block_n_strains(self) -> pd.Series:
        """Returns a pandas series with columns=blocks, values=#paths in which the block is present"""
        return self.df.groupby("block_id")["path_id"].nunique()

    def to_blockstats_df(self) -> pd.DataFrame:
        """Returns a dataframe containing statistics about blocks distribution.
        - count: n. times the block occurs
        - n. strains: number of strains in which the block is observed
        - duplicated: whether the block is duplicated in at least one strain
        - core: whether a block occurrs exactly once per strain
        """
        block_counts = self.block_counts()
        block_n_strains = self.block_n_strains()

        df = pd.DataFrame({"count": block_counts, "n_strains": block_n_strains})
        df["duplicated"] = df["count"] > df["n_strains"]
        df["core"] = df["n_strains"] == len(self.df["path_id"].unique())
        df["core"] &= ~df["duplicated"]
        return df

    def node_to_block(self, node_id: int) -> tuple[int, bool]:
        """Returns the node's block id and strandedness"""
        return self.df.loc[str(node_id), ["block_id", "strand"]]

    def nodes_to_blocks(self, node_ids: list[int]) -> tuple[list[int], list[bool]]:
        """Returns the block id and strandedness of a list of nodes"""
        N = np.array(node_ids, dtype=str)
        return self.df.loc[N, ["block_id", "strand"]].values.T
