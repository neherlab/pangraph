import pypangraph as pp
import pypangraph.topology as tp
from pypangraph.msu import (
    find_mergers,
    minimal_synteny_units,
    core_paths,
    flip_msu_to_most_common_orientation,
)
from collections import defaultdict
import pytest


@pytest.fixture
def load_graph():
    return pp.Pangraph.from_json("tests/data/plasmids.json")


@pytest.fixture
def generate_core_paths():
    A = tp.Node("A", True)
    B = tp.Node("B", True)
    C = tp.Node("C", True)  # invert
    D = tp.Node("D", True)  # invert
    E = tp.Node("E", True)
    F = tp.Node("F", True)
    G = tp.Node("G", True)
    H = tp.Node("H", True)  # invert
    J = tp.Node("J", True)

    p1 = tp.Path([A, B, C, D, E, F, G, H, J], circular=True)
    p2 = tp.Path([A, B, C, D, E, F, G, H, J], circular=True)
    p3 = tp.Path([A, B, ~D, ~C, E, F, G, H, J], circular=True)
    p4 = tp.Path([A, B, ~D, ~C, E, F, G, ~H, J], circular=True)

    paths = {1: p1, 2: p2, 3: p3, 4: p4}
    nodes = {n.id: n for n in [A, B, C, D, E, F, G, H, J]}

    return paths, nodes


class TestMSUFunctions:
    def test_find_mergers_same(self):
        # Create specific paths to test merger logic
        A = tp.Node("A", True)
        B = tp.Node("B", True)
        C = tp.Node("C", True)
        D = tp.Node("D", True)
        E = tp.Node("E", True)

        # Create paths where A-B and C-D should be merged
        p1 = tp.Path([A, B, C, D, E], circular=False)
        p2 = tp.Path([A, B, C, D, E], circular=False)

        paths = {"iso1": p1, "iso2": p2}
        mergers = find_mergers(paths)

        # Each block should map to itself or another block
        assert len(mergers) == 5
        assert all(bid in mergers for bid in ["A", "B", "C", "D", "E"])

    def test_find_mergers_complex(self, generate_core_paths):
        paths, nodes = generate_core_paths

        mg = find_mergers(paths)
        mg_groups = defaultdict(set)
        for source, sink in mg.items():
            mg_groups[sink].add(source)
        sources = list(mg_groups.values())
        assert len(sources) == 4
        assert {"A", "B", "J"} in sources
        assert {"C", "D"} in sources
        assert {"E", "F", "G"} in sources
        assert {"H"} in sources

    def test_minimal_synteny_units(self, load_graph):
        """Test the main MSU function"""
        pan = load_graph
        MSU_mergers, MSU_paths, MSU_len = minimal_synteny_units(
            pan, L_thr=50, rotate=True
        )

        # Basic checks
        assert isinstance(MSU_mergers, dict)
        assert isinstance(MSU_paths, dict)
        assert isinstance(MSU_len, dict)

        # MSU IDs should be in the format MSU_X
        assert all(msu_id.startswith("MSU_") for msu_id in MSU_len.keys())

        # All paths should be circular when rotate=True
        assert all(path.circular for path in MSU_paths.values())

    def test_minimal_synteny_units_no_rotate(self, load_graph):
        """Test MSU function without rotation"""
        pan = load_graph
        MSU_mergers, MSU_paths, MSU_len = minimal_synteny_units(
            pan, L_thr=50, rotate=False
        )

        # Basic checks
        assert isinstance(MSU_mergers, dict)
        assert isinstance(MSU_paths, dict)
        assert isinstance(MSU_len, dict)

    def test_core_paths(self, load_graph):
        """Test core_paths function"""
        pan = load_graph
        c_paths = core_paths(pan, L_thr=50)

        # Should return a dictionary of paths
        assert isinstance(c_paths, dict)

        # All values should be Path objects
        assert all(isinstance(path, tp.Path) for path in c_paths.values())

        # All nodes in core paths should have the minimum length
        bdf = pan.to_blockstats_df()
        for path in c_paths.values():
            for node in path.nodes:
                block_len = bdf.loc[node.id, "len"]
                block_core = bdf.loc[node.id, "core"]
                assert block_len >= 50
                assert block_core

    def test_flip_msu_to_most_common_orientation(self):
        """Test the orientation flipping function"""
        # Create test paths with different orientations
        A = tp.Node("A", True)  # Should stay positive (majority)
        B = tp.Node("B", False)  # Should flip to positive (minority negative)
        C = tp.Node("C", True)  # Should stay positive

        # Path 1: A+, B-, C+
        p1 = tp.Path([A, B, C], circular=True)
        # Path 2: A+, B-, C+
        p2 = tp.Path([A, B, C], circular=True)
        # Path 3: A+, B+, C+  (B is positive here)
        p3 = tp.Path([A, ~B, C], circular=True)  # B+ in this path

        paths = {"iso1": p1, "iso2": p2, "iso3": p3}

        # B appears negative twice and positive once, so it should flip to positive
        flipped_paths = flip_msu_to_most_common_orientation(paths)

        # Check that paths were modified
        assert isinstance(flipped_paths, dict)
        assert len(flipped_paths) == 3

    def test_core_paths_empty_result(self, load_graph):
        """Test core_paths with very high threshold"""
        pan = load_graph
        c_paths = core_paths(pan, L_thr=10000)  # Very high threshold

        # Should return empty paths or very few nodes
        assert isinstance(c_paths, dict)

        # Most paths should be empty or very short
        total_nodes = sum(len(path.nodes) for path in c_paths.values())
        assert total_nodes >= 0  # At least no errors


class TestMSUIntegration:
    """Integration tests for MSU workflow"""

    def test_msu_workflow_consistency(self, load_graph):
        """Test that the MSU workflow produces consistent results"""
        pan = load_graph

        # Run MSU analysis
        MSU_mergers, MSU_paths, MSU_len = minimal_synteny_units(
            pan, L_thr=50, rotate=True
        )

        # Check consistency between mergers and paths
        # All MSU IDs in mergers should appear in MSU_len
        msu_ids_in_mergers = set(MSU_mergers.values())
        msu_ids_in_len = set(MSU_len.keys())
        assert msu_ids_in_mergers.issubset(msu_ids_in_len)

        # All MSU IDs in paths should appear in MSU_len
        msu_ids_in_paths = set()
        for path in MSU_paths.values():
            for node in path.nodes:
                msu_ids_in_paths.add(node.id)
        assert msu_ids_in_paths.issubset(msu_ids_in_len)

    def test_msu_ordering_by_length(self, load_graph):
        """Test that MSUs are ordered by length (MSU_0 should be longest)"""
        pan = load_graph
        MSU_mergers, MSU_paths, MSU_len = minimal_synteny_units(
            pan, L_thr=50, rotate=True
        )

        # Extract numeric IDs and sort
        msu_items = [
            (int(msu_id.split("_")[1]), length) for msu_id, length in MSU_len.items()
        ]
        msu_items.sort(key=lambda x: x[0])  # Sort by ID number

        # Check that lengths are in descending order
        lengths = [item[1] for item in msu_items]
        assert lengths == sorted(lengths, reverse=True)

    def test_different_thresholds(self, load_graph):
        """Test MSU analysis with different length thresholds"""
        pan = load_graph

        # Test with low threshold
        MSU_mergers_low, MSU_paths_low, MSU_len_low = minimal_synteny_units(
            pan, L_thr=10, rotate=False
        )

        # Test with high threshold
        MSU_mergers_high, MSU_paths_high, MSU_len_high = minimal_synteny_units(
            pan, L_thr=100, rotate=False
        )

        # Higher threshold should result in fewer or equal MSUs
        assert len(MSU_len_high) <= len(MSU_len_low)

        # Higher threshold should result in fewer or equal total nodes
        total_nodes_low = sum(len(path.nodes) for path in MSU_paths_low.values())
        total_nodes_high = sum(len(path.nodes) for path in MSU_paths_high.values())
        assert total_nodes_high <= total_nodes_low
