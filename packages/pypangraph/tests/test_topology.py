import pypangraph.topology as tp
from collections import defaultdict
import pytest


class TestNode:
    def test_eq(self):
        n1 = tp.Node(1, True)
        n2 = tp.Node(1, True)
        n3 = tp.Node(1, False)
        n4 = tp.Node(2, True)
        assert n1 == n2
        assert n1 != n3
        assert n1 != n4

    def test_invert(self):
        n = tp.Node(1, True)
        assert n.invert() == tp.Node(1, False)

    def test_hash(self):
        n1 = tp.Node(1, True)
        n2 = tp.Node(1, True)
        n3 = tp.Node(1, False)
        assert hash(n1) == hash(n2)
        assert hash(n1) != hash(n3)

    def test_repr(self):
        n1 = tp.Node(1, True)
        n2 = tp.Node(1, False)
        assert str(n1) == "[1|+]"
        assert str(n2) == "[1|-]"

    def test_to_str_id(self):
        n1 = tp.Node(1, True)
        n2 = tp.Node(1, False)
        assert n1.to_str_id() == "1_f"
        assert n2.to_str_id() == "1_r"

    def test_from_str_id(self):
        n1 = tp.Node.from_str_id("1_f")
        n2 = tp.Node.from_str_id("1_r")
        assert n1 == tp.Node(1, True)
        assert n2 == tp.Node(1, False)

    def test_invert_operator(self):
        n = tp.Node(1, True)
        assert ~n == tp.Node(1, False)


class TestPath:
    def test_eq(self):
        n1 = tp.Node(1, True)
        n2 = tp.Node(2, True)
        n3 = tp.Node(3, True)
        p1 = tp.Path([n1, n2, n3], circular=True)
        p2 = tp.Path([n1, n2, n3], circular=True)
        p3 = tp.Path([n1, ~n2, n3], circular=True)
        p4 = tp.Path([n1, n2], circular=True)
        assert p1 == p2
        assert p1 != p3
        assert p1 != p4

    def test_invert(self):
        n1 = tp.Node(1, True)
        n2 = tp.Node(2, True)
        n3 = tp.Node(3, True)
        p = tp.Path([n1, n2, n3], circular=True)
        assert ~p == tp.Path([~n3, ~n2, ~n1], circular=True)

    def test_rotate_to(self):
        n1 = tp.Node(1, True)
        n2 = tp.Node(2, True)
        n3 = tp.Node(3, True)
        n4 = tp.Node(4, True)
        p = tp.Path([n1, n2, n3, n4], circular=True)
        assert p.rotate_to(2, True) == tp.Path([n2, n3, n4, n1], circular=True)
        assert p.rotate_to(4, True) == tp.Path([n4, n1, n2, n3], circular=True)
        assert p.rotate_to(2, False) == tp.Path([~n2, ~n1, ~n4, ~n3], circular=True)
        assert p.rotate_to(4, False) == tp.Path([~n4, ~n3, ~n2, ~n1], circular=True)

    def test_rotate_to_node(self):
        n1 = tp.Node(1, True)
        n2 = tp.Node(2, True)
        n3 = tp.Node(3, True)
        p = tp.Path([n1, n2, n3], circular=True)
        assert p.rotate_to_node(n2) == tp.Path([n2, n3, n1], circular=True)

    def test_rename_bids(self):
        n1 = tp.Node(1, True)
        n2 = tp.Node(2, False)
        n3 = tp.Node(3, True)
        p = tp.Path([n1, n2, n3], circular=True)
        assert p.rename_bids({1: 10, 2: 20, 3: 30}) == tp.Path(
            [tp.Node(10, True), tp.Node(20, False), tp.Node(30, True)],
            circular=True,
        )

    def test_add_left(self):
        n1 = tp.Node(1, True)
        n2 = tp.Node(2, True)
        p = tp.Path([n1], circular=False)
        p.add_left(n2)
        assert p.nodes == [n2, n1]

    def test_add_right(self):
        n1 = tp.Node(1, True)
        n2 = tp.Node(2, True)
        p = tp.Path([n1], circular=False)
        p.add_right(n2)
        assert p.nodes == [n1, n2]

    def test_len(self):
        n1 = tp.Node(1, True)
        n2 = tp.Node(2, True)
        p = tp.Path([n1, n2], circular=False)
        assert len(p) == 2

    def test_hash(self):
        n1 = tp.Node(1, True)
        n2 = tp.Node(2, True)
        p1 = tp.Path([n1, n2], circular=True)
        p2 = tp.Path([n1, n2], circular=True)
        assert hash(p1) == hash(p2)

    def test_repr(self):
        n1 = tp.Node(1, True)
        n2 = tp.Node(2, True)
        p = tp.Path([n1, n2], circular=False)
        assert str(p) == "[1|+]_[2|+]"

    def test_to_list(self):
        n1 = tp.Node(1, True)
        n2 = tp.Node(2, False)
        p = tp.Path([n1, n2], circular=False)
        assert p.to_list() == ["1_f", "2_r"]

    def test_from_list(self):
        p = tp.Path.from_list(["1_f", "2_r"], circular=False)
        expected = tp.Path([tp.Node(1, True), tp.Node(2, False)], circular=False)
        assert p.nodes[0] == expected.nodes[0]
        assert p.nodes[1] == expected.nodes[1]

    def test_rotate_to_assertion_errors(self):
        n1 = tp.Node(1, True)
        n2 = tp.Node(2, True)
        p = tp.Path([n1, n2], circular=False)

        # Should raise assertion error for non-circular path
        with pytest.raises(AssertionError, match="Path is not circular"):
            p.rotate_to(1, True)

        # Should raise assertion error for block not in path
        p_circular = tp.Path([n1, n2], circular=True)
        with pytest.raises(AssertionError, match="Block not in path"):
            p_circular.rotate_to(3, True)


class TestEdge:
    def test_invert(self):
        n1 = tp.Node(1, True)
        n2 = tp.Node(2, True)
        e = tp.Edge(n1, n2)
        assert e.invert() == tp.Edge(n2.invert(), n1.invert())

    def test_eq(self):
        n1 = tp.Node(1, True)
        n2 = tp.Node(2, True)
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
        n1 = tp.Node(1, True)
        n2 = tp.Node(2, True)
        e1 = tp.Edge(n1, n2)
        e2 = tp.Edge(n1, n2)  # Same edge
        e3 = tp.Edge(n1, ~n2)  # Different edge
        assert hash(e1) == hash(e2)  # Same edges should have same hash
        assert hash(e1) != hash(e3)  # Different edges should have different hash

    def test_repr(self):
        n1 = tp.Node(1, True)
        n2 = tp.Node(2, True)
        e = tp.Edge(n1, n2)
        assert str(e) == "[1|+] <--> [2|+]"

    def test_invert_operator(self):
        n1 = tp.Node(1, True)
        n2 = tp.Node(2, True)
        e = tp.Edge(n1, n2)
        assert ~e == tp.Edge(n2.invert(), n1.invert())

    def test_to_str_id(self):
        n1 = tp.Node(1, True)
        n2 = tp.Node(2, True)
        e = tp.Edge(n1, n2)
        str_id = e.to_str_id()
        assert str_id in ["1_f__2_f", "2_r__1_r"]  # Should be consistent ordering

    def test_from_str_id(self):
        e = tp.Edge.from_str_id("1_f__2_f")
        expected = tp.Edge(tp.Node(1, True), tp.Node(2, True))
        assert e.left == expected.left
        assert e.right == expected.right


class TestTopologyFunctions:
    @pytest.fixture
    def sample_paths(self):
        """Create sample paths for testing"""
        A = tp.Node(1, True)
        B = tp.Node(2, True)
        C = tp.Node(3, True)
        D = tp.Node(4, True)
        E = tp.Node(5, True)

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
            return bid in [1, 2, 3]

        filtered = tp.filter_paths(sample_paths, keep_A_B_C)

        # Should keep 1, 2, 3 but filter out 4, 5
        assert len(filtered["iso1"].nodes) == 3
        assert len(filtered["iso2"].nodes) == 3
        assert len(filtered["iso3"].nodes) == 2  # 1, 2 (4 filtered out)
        assert (
            len(filtered["iso4"].nodes) == 0
        )  # 4, 5 filtered out    def test_path_categories(self, sample_paths):
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

        # 1->2 edge should appear in multiple paths
        A_B_edge = tp.Edge(tp.Node(1, True), tp.Node(2, True))
        assert A_B_edge in edge_count
        assert edge_count[A_B_edge] == 3

    def test_path_block_count(self, sample_paths):
        block_count = tp.path_block_count(sample_paths)

        # Check that blocks are counted correctly
        assert block_count[1] == 3  # Appears in iso1, iso2, iso3
        assert block_count[2] == 3  # Appears in iso1, iso2, iso3
        assert block_count[3] == 2  # Appears in iso1, iso2
        assert block_count[4] == 2  # Appears in iso3, iso4
        assert block_count[5] == 1  # Appears in iso4
