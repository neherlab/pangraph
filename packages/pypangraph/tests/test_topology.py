import pypangraph.topology as tp
from collections import defaultdict
import pytest


class TestNode:
    def test_eq(self):
        n1 = tp.Node("A", True)
        n2 = tp.Node("A", True)
        n3 = tp.Node("A", False)
        n4 = tp.Node("B", True)
        assert n1 == n2
        assert n1 != n3
        assert n1 != n4

    def test_invert(self):
        n = tp.Node("A", True)
        assert n.invert() == tp.Node("A", False)

    def test_hash(self):
        n1 = tp.Node("A", True)
        n2 = tp.Node("A", True)
        n3 = tp.Node("A", False)
        assert hash(n1) == hash(n2)
        assert hash(n1) != hash(n3)

    def test_repr(self):
        n1 = tp.Node("A", True)
        n2 = tp.Node("A", False)
        assert str(n1) == "[A|+]"
        assert str(n2) == "[A|-]"

    def test_to_str_id(self):
        n1 = tp.Node("A", True)
        n2 = tp.Node("A", False)
        assert n1.to_str_id() == "A_f"
        assert n2.to_str_id() == "A_r"

    def test_from_str_id(self):
        n1 = tp.Node.from_str_id("A_f")
        n2 = tp.Node.from_str_id("A_r")
        assert n1 == tp.Node("A", True)
        assert n2 == tp.Node("A", False)

    def test_invert_operator(self):
        n = tp.Node("A", True)
        assert ~n == tp.Node("A", False)


class TestPath:
    def test_eq(self):
        n1 = tp.Node("A", True)
        n2 = tp.Node("B", True)
        n3 = tp.Node("C", True)
        p1 = tp.Path([n1, n2, n3], circular=True)
        p2 = tp.Path([n1, n2, n3], circular=True)
        p3 = tp.Path([n1, ~n2, n3], circular=True)
        p4 = tp.Path([n1, n2], circular=True)
        assert p1 == p2
        assert p1 != p3
        assert p1 != p4

    def test_invert(self):
        n1 = tp.Node("A", True)
        n2 = tp.Node("B", True)
        n3 = tp.Node("C", True)
        p = tp.Path([n1, n2, n3], circular=True)
        assert ~p == tp.Path([~n3, ~n2, ~n1], circular=True)

    def test_rotate_to(self):
        n1 = tp.Node("A", True)
        n2 = tp.Node("B", True)
        n3 = tp.Node("C", True)
        n4 = tp.Node("D", True)
        p = tp.Path([n1, n2, n3, n4], circular=True)
        assert p.rotate_to("B", True) == tp.Path([n2, n3, n4, n1], circular=True)
        assert p.rotate_to("D", True) == tp.Path([n4, n1, n2, n3], circular=True)
        assert p.rotate_to("B", False) == tp.Path([~n2, ~n1, ~n4, ~n3], circular=True)
        assert p.rotate_to("D", False) == tp.Path([~n4, ~n3, ~n2, ~n1], circular=True)

    def test_rotate_to_node(self):
        n1 = tp.Node("A", True)
        n2 = tp.Node("B", True)
        n3 = tp.Node("C", True)
        p = tp.Path([n1, n2, n3], circular=True)
        assert p.rotate_to_node(n2) == tp.Path([n2, n3, n1], circular=True)

    def test_rename_bids(self):
        n1 = tp.Node("A", True)
        n2 = tp.Node("B", False)
        n3 = tp.Node("C", True)
        p = tp.Path([n1, n2, n3], circular=True)
        assert p.rename_bids({"A": "X", "B": "Y", "C": "Z"}) == tp.Path(
            [tp.Node("X", True), tp.Node("Y", False), tp.Node("Z", True)],
            circular=True,
        )

    def test_add_left(self):
        n1 = tp.Node("A", True)
        n2 = tp.Node("B", True)
        p = tp.Path([n1], circular=False)
        p.add_left(n2)
        assert p.nodes == [n2, n1]

    def test_add_right(self):
        n1 = tp.Node("A", True)
        n2 = tp.Node("B", True)
        p = tp.Path([n1], circular=False)
        p.add_right(n2)
        assert p.nodes == [n1, n2]

    def test_len(self):
        n1 = tp.Node("A", True)
        n2 = tp.Node("B", True)
        p = tp.Path([n1, n2], circular=False)
        assert len(p) == 2

    def test_hash(self):
        n1 = tp.Node("A", True)
        n2 = tp.Node("B", True)
        p1 = tp.Path([n1, n2], circular=True)
        p2 = tp.Path([n1, n2], circular=True)
        assert hash(p1) == hash(p2)

    def test_repr(self):
        n1 = tp.Node("A", True)
        n2 = tp.Node("B", True)
        p = tp.Path([n1, n2], circular=False)
        assert str(p) == "[A|+]_[B|+]"

    def test_to_list(self):
        n1 = tp.Node("A", True)
        n2 = tp.Node("B", False)
        p = tp.Path([n1, n2], circular=False)
        assert p.to_list() == ["A_f", "B_r"]

    def test_from_list(self):
        p = tp.Path.from_list(["A_f", "B_r"], circular=False)
        expected = tp.Path([tp.Node("A", True), tp.Node("B", False)], circular=False)
        assert p.nodes[0] == expected.nodes[0]
        assert p.nodes[1] == expected.nodes[1]

    def test_rotate_to_assertion_errors(self):
        n1 = tp.Node("A", True)
        n2 = tp.Node("B", True)
        p = tp.Path([n1, n2], circular=False)

        # Should raise assertion error for non-circular path
        with pytest.raises(AssertionError, match="Path is not circular"):
            p.rotate_to("A", True)

        # Should raise assertion error for block not in path
        p_circular = tp.Path([n1, n2], circular=True)
        with pytest.raises(AssertionError, match="Block not in path"):
            p_circular.rotate_to("C", True)


class TestEdge:
    def test_invert(self):
        n1 = tp.Node("A", True)
        n2 = tp.Node("B", True)
        e = tp.Edge(n1, n2)
        assert e.invert() == tp.Edge(n2.invert(), n1.invert())

    def test_eq(self):
        n1 = tp.Node("A", True)
        n2 = tp.Node("B", True)
        e1 = tp.Edge(n1, n2)
        e2 = tp.Edge(n1, n2)
        e3 = tp.Edge(~n1, ~n2)
        e4 = tp.Edge(n1, ~n2)
        e5 = tp.Edge(~n2, ~n1)
        assert e1 == e2
        assert e1 != e3
        assert e1 != e4
        assert e1 == e5

    def test_hash(self):
        n1 = tp.Node("A", True)
        n2 = tp.Node("B", True)
        e1 = tp.Edge(n1, n2)
        e2 = tp.Edge(n1, n2)  # Same edge
        e3 = tp.Edge(n1, ~n2)  # Different edge
        assert hash(e1) == hash(e2)  # Same edges should have same hash
        assert hash(e1) != hash(e3)  # Different edges should have different hash

    def test_repr(self):
        n1 = tp.Node("A", True)
        n2 = tp.Node("B", True)
        e = tp.Edge(n1, n2)
        assert str(e) == "[A|+] <--> [B|+]"

    def test_invert_operator(self):
        n1 = tp.Node("A", True)
        n2 = tp.Node("B", True)
        e = tp.Edge(n1, n2)
        assert ~e == tp.Edge(n2.invert(), n1.invert())

    def test_to_str_id(self):
        n1 = tp.Node("A", True)
        n2 = tp.Node("B", True)
        e = tp.Edge(n1, n2)
        str_id = e.to_str_id()
        assert str_id in ["A_f__B_f", "B_r__A_r"]  # Should be consistent ordering

    def test_from_str_id(self):
        e = tp.Edge.from_str_id("A_f__B_f")
        expected = tp.Edge(tp.Node("A", True), tp.Node("B", True))
        assert e.left == expected.left
        assert e.right == expected.right


class TestTopologyFunctions:
    @pytest.fixture
    def sample_paths(self):
        """Create sample paths for testing"""
        A = tp.Node("A", True)
        B = tp.Node("B", True)
        C = tp.Node("C", True)
        D = tp.Node("D", True)
        E = tp.Node("E", True)

        p1 = tp.Path([A, B, C], circular=True)
        p2 = tp.Path([A, B, C], circular=True)
        p3 = tp.Path([A, B, D], circular=True)
        p4 = tp.Path([D, E], circular=False)

        return {
            "iso1": p1,
            "iso2": p2,
            "iso3": p3,
            "iso4": p4,
        }

    def test_filter_paths(self, sample_paths):
        def keep_A_B_C(bid):
            return bid in ["A", "B", "C"]

        filtered = tp.filter_paths(sample_paths, keep_A_B_C)

        # Should keep A, B, C but filter out D, E
        assert len(filtered["iso1"].nodes) == 3
        assert len(filtered["iso2"].nodes) == 3
        assert len(filtered["iso3"].nodes) == 2  # A, B (D filtered out)
        assert len(filtered["iso4"].nodes) == 0  # D, E filtered out

    def test_path_categories(self, sample_paths):
        categories = tp.path_categories(sample_paths)

        # Should be sorted by count (descending)
        assert len(categories) == 3  # 3 unique non-empty paths

        # First category should be the most common (appears twice)
        count, nodes, isolates = categories[0]
        assert count == 2
        assert len(isolates) == 2

    def test_path_edge_count(self, sample_paths):
        edge_count = tp.path_edge_count(sample_paths)

        # Check that edges are counted correctly
        assert len(edge_count) > 0

        # A->B edge should appear in multiple paths
        A_B_edge = tp.Edge(tp.Node("A", True), tp.Node("B", True))
        assert A_B_edge in edge_count
        assert edge_count[A_B_edge] == 3

    def test_path_block_count(self, sample_paths):
        block_count = tp.path_block_count(sample_paths)

        # Check that blocks are counted correctly
        assert block_count["A"] == 3  # Appears in iso1, iso2, iso3
        assert block_count["B"] == 3  # Appears in iso1, iso2, iso3
        assert block_count["C"] == 2  # Appears in iso1, iso2
        assert block_count["D"] == 2  # Appears in iso3, iso4
        assert block_count["E"] == 1  # Appears in iso4
