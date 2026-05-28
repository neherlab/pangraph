import pandas as pd
from Bio.SeqRecord import SeqRecord

from ..topology_utils import Walk
from .junction import JunctionNode, Junction, path_junction_split
from .dataframe import junctions_dataframe
from .positions import junction_positions
from .sequences import junction_sequences
from .stats import junction_stats


class BackboneJunctions:
    """Backbone junction analysis for a pangenome graph.

    Splits each path into junctions at backbone (core + above-threshold length)
    block boundaries. Results are lazily computed and cached.

    Args:
        pan: A Pangraph object.
        L_thr: Minimum block length to be considered backbone (default 500 bp).
    """

    def __init__(self, pan, L_thr: int = 500):
        self.pan = pan
        self.L_thr = L_thr
        self._bdf = pan.to_blockstats_df()
        self._junctions = None  # dict[iso, list[Junction]]
        self._edge_map = None  # dict[edge_str, dict[iso, Junction]]

    def _is_backbone(self, bid: str) -> bool:
        """Determine if a block ID corresponds to a backbone block."""
        return self._bdf.loc[bid, "core"] and self._bdf.loc[bid, "len"] >= self.L_thr

    def _ensure_split(self):
        """Ensure that paths have been split into junctions and edge maps built."""
        if self._junctions is not None:
            return
        self._junctions = {}
        self._edge_map = {}
        for name, path in self.pan.paths.items():
            # create a Walk from a Path object, including accessory blocks
            nodes = []
            for nid in path.nodes:
                row = self.pan.nodes.df.loc[str(nid)]
                nodes.append(JunctionNode(row["block_id"], row["strand"], nid))
            tu_path = Walk(nodes, path.circular)

            juncs = path_junction_split(tu_path, self._is_backbone)
            self._junctions[name] = juncs
            for j in juncs:
                edge = j.flanking_edge()
                if edge is None:
                    continue
                edge_str = edge.to_str_id()
                if edge_str not in self._edge_map:
                    self._edge_map[edge_str] = {}
                self._edge_map[edge_str][name] = j

    def junctions_for(self, isolate: str) -> list[Junction]:
        """Return all junctions for a given isolate."""
        self._ensure_split()
        return self._junctions[isolate]

    def junction_for(self, isolate: str, edge_str: str) -> Junction:
        """Return the junction for a given isolate and edge."""
        self._ensure_split()
        return self._edge_map[edge_str][isolate]

    def edges(self) -> list[str]:
        """Return list of all edge string IDs."""
        self._ensure_split()
        return list(self._edge_map.keys())

    def stats(self) -> pd.DataFrame:
        """Compute per-edge junction statistics.

        Returns:
            DataFrame with edge string IDs as index and columns:
            frequency, n_categories, majority_category_freq, is_transitive,
            is_singleton, left_core_length, right_core_length, accessory_length.
            Sorted by frequency descending.
        """
        self._ensure_split()
        return junction_stats(self._edge_map, self._bdf)

    def dataframe(self) -> tuple[pd.DataFrame, pd.DataFrame]:
        """Build a pivot table of junction lengths per isolate/edge.

        Returns:
            A tuple (jdf, stats_df) where:
            - jdf: DataFrame with isolates as rows, edges as columns, accessory
              lengths as values. NaN for absent junctions. Columns sorted by
              frequency descending.
            - stats_df: Per-edge statistics DataFrame (see stats() for columns).
        """
        self._ensure_split()
        stats_df = self.stats()
        jdf = junctions_dataframe(self._edge_map, self._bdf)
        if not jdf.empty:
            # Sort columns by frequency (same order as stats_df index)
            jdf = jdf[stats_df.index]
        return jdf, stats_df

    def positions(self) -> pd.DataFrame:
        """Find genomic positions of flanking core blocks for each junction.

        Returns:
            A DataFrame with MultiIndex (iso, edge) and columns:
            - left_start, left_end: genomic position of the left flanking block
            - right_start, right_end: genomic position of the right flanking block
            - strand: True if canonical orientation, False if inverted
        """
        self._ensure_split()
        return junction_positions(self._edge_map, self.pan)

    def sequences(self, edge_str: str) -> list[SeqRecord]:
        """Extract co-oriented sequences spanning a junction.

        For each isolate with the given junction, returns a SeqRecord spanning
        from the start of the left core block to the end of the right one.
        All sequences are co-oriented: core blocks are in the same orientation
        across isolates, with accessory sequence in between. Individual blocks
        are reverse-complemented as needed based on their strand.

        Args:
            edge_str: The canonical edge string ID (e.g. "100_f__200_f").

        Returns:
            A list of SeqRecord objects, one per isolate. The record id is
            the isolate name, and the description contains the edge string ID.
        """
        self._ensure_split()
        return junction_sequences(self._edge_map, self.pan, edge_str)
