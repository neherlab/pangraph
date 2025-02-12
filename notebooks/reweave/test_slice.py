from slice_utils import *
import unittest


class TestSlice(unittest.TestCase):

    def generate_example(self):
        seq = "ACTGGATATCCGATATTCGAG"
        s = lambda pos, base: Substitution(pos, base)
        d = lambda pos, length: Deletion(pos, length)
        i = lambda pos, ins: Insertion(pos, ins)
        ed = Edit(
            subs=[
                s(2, "C"),
                s(5, "C"),
                s(6, "G"),
                s(7, "C"),
                s(13, "G"),
                s(14, "T"),
                s(18, "C"),
                s(20, "A"),
            ],
            dels=[
                d(0, 2),
                d(4, 3),
                d(9, 2),
                d(13, 4),
                d(18, 3),
            ],
            ins=[
                i(2, "CC"),
                i(5, "A"),
                i(6, "TTT"),
                i(10, "C"),
                i(13, "T"),
                i(14, "GG"),
                i(17, "A"),
                i(21, "A"),
            ],
        )
        return seq, ed

    def test_slice_substitutions(self):
        seq, ed = self.generate_example()
        i = Interval(start=6, end=14, aligned=True, new_block_id=0)
        new_ed = slice_substitutions(i, ed.subs)
        self.assertEqual(
            new_ed,
            [
                Substitution(0, "G"),
                Substitution(1, "C"),
                Substitution(7, "G"),
            ],
        )
        i = Interval(start=15, end=21, aligned=True, new_block_id=0)
        new_ed = slice_substitutions(i, ed.subs)
        self.assertEqual(
            new_ed,
            [
                Substitution(3, "C"),
                Substitution(5, "A"),
            ],
        )

    def test_slice_deletions(self):
        seq, ed = self.generate_example()
        i = Interval(start=6, end=14, aligned=True, new_block_id=0)
        new_ed = slice_deletions(i, ed.dels)
        self.assertEqual(
            new_ed,
            [
                Deletion(0, 1),
                Deletion(3, 2),
                Deletion(7, 1),
            ],
        )
        i = Interval(start=15, end=21, aligned=True, new_block_id=0)
        new_ed = slice_deletions(i, ed.dels)
        self.assertEqual(
            new_ed,
            [
                Deletion(0, 2),
                Deletion(3, 3),
            ],
        )

    def test_slice_insertions(self):
        seq, ed = self.generate_example()
        i = Interval(start=6, end=14, aligned=True, new_block_id=0)
        new_ed = slice_insertions(i, ed.ins, len(seq))
        self.assertEqual(
            new_ed,
            [
                Insertion(0, "TTT"),
                Insertion(4, "C"),
                Insertion(7, "T"),
            ],
        )
        i = Interval(start=15, end=21, aligned=True, new_block_id=0)
        new_ed = slice_insertions(i, ed.ins, len(seq))
        self.assertEqual(
            new_ed,
            [
                Insertion(2, "A"),
                Insertion(6, "A"),
            ],
        )

    def test_interval_node_coords(self):
        seq, ed = self.generate_example()

        i = Interval(start=6, end=14, aligned=True, new_block_id=0)
        new_pos = interval_node_coords(i, ed, len(seq))
        self.assertEqual(new_pos, (5, 14))

        i = Interval(start=15, end=21, aligned=True, new_block_id=0)
        new_pos = interval_node_coords(i, ed, len(seq))
        self.assertEqual(new_pos, (16, 19))


class TestPosition(unittest.TestCase):

    def test_new_position(self):
        path_L = 100

        strandedness = True
        node_coords = (10, 20)
        old_position = (10, 40)
        new_pos = new_position(old_position, node_coords, path_L, strandedness)
        self.assertEqual(new_pos, (20, 30))

        old_position = (95, 20)
        new_pos = new_position(old_position, node_coords, path_L, strandedness)
        self.assertEqual(new_pos, (5, 15))

        strandedness = False
        old_position = (10, 50)
        new_pos = new_position(old_position, node_coords, path_L, strandedness)
        self.assertEqual(new_pos, (30, 40))

        old_position = (40, 5)
        new_pos = new_position(old_position, node_coords, path_L, strandedness)
        self.assertEqual(new_pos, (85, 95))

    def test_node_coords(self):
        i = Interval(start=10, end=20, aligned=True, new_block_id=0)
        ed = Edit(
            subs=[Substitution(2, "G"), Substitution(13, "T"), Substitution(24, "T")],
            dels=[Deletion(18, 3)],
            ins=[Insertion(7, "A"), Insertion(10, "AAAA"), Insertion(20, "TTTTTTTT")],
        )
        block_L = 100
        new_pos = interval_node_coords(i, ed, block_L)
        self.assertEqual(new_pos, (11, 23))


class TestBlockSlice(unittest.TestCase):

    def generate_example(self):
        seq = "ACTTGATCCTTATATTTATCCGATCAT"
        bid = 1
        s = lambda pos, base: Substitution(pos, base)
        d = lambda pos, length: Deletion(pos, length)
        i = lambda pos, ins: Insertion(pos, ins)
        ed1 = Edit(
            subs=[s(2, "G"), s(13, "T"), s(24, "T")],
            dels=[d(18, 3)],
            ins=[i(7, "A"), i(10, "A")],
        )
        ed2 = Edit(
            subs=[s(4, "T"), s(19, "G"), s(20, "G")],
            dels=[d(6, 2), d(13, 2)],
            ins=[i(17, "T"), i(25, "A")],
        )
        ed3 = Edit(
            subs=[],
            dels=[d(2, 4), d(9, 3), d(24, 2)],
            ins=[i(20, "T")],
        )
        n1 = Node(1, bid, path_id=1, position=(100, 125), strandedness=True)
        n2 = Node(2, bid, path_id=2, position=(1000, 1025), strandedness=False)
        n3 = Node(3, bid, path_id=3, position=(90, 9), strandedness=False)
        p1 = Path(1, [1, 4], L=2000)
        p2 = Path(2, [2, 5], L=2000)
        p3 = Path(3, [3, 6], L=100)
        b1 = Block(id=bid, consensus=seq, alignment={1: ed1, 2: ed2, 3: ed3})
        G = Pangraph(
            paths={1: p1, 2: p2, 3: p3}, blocks={bid: b1}, nodes={1: n1, 2: n2, 3: n3}
        )
        return b1, G

    def test_new_strandedness(self):
        ns = lambda d, o, s: new_strandedness(
            old_strandedness=s, orientation=o, is_anchor=d
        )
        self.assertEqual(ns(True, True, True), True)
        self.assertEqual(ns(True, True, False), False)
        self.assertEqual(ns(True, False, True), True)
        self.assertEqual(ns(True, False, False), False)
        self.assertEqual(ns(False, True, True), True)
        self.assertEqual(ns(False, True, False), False)
        self.assertEqual(ns(False, False, True), False)
        self.assertEqual(ns(False, False, False), True)

    def test_block_slice_fwd_anchor(self):
        b, G = self.generate_example()
        new_bid = 42
        i = Interval(
            start=10,
            end=20,
            aligned=True,
            new_block_id=new_bid,
            orientation=True,
            is_anchor=True,
        )

        new_b, new_nodes = block_slice(b, i, G)

        # block consensus
        self.assertEqual(new_b.consensus, "TATATTTATC")

        # new nodes and genome position
        pos1 = (111, 120)
        nn1 = Node(None, block_id=new_bid, path_id=1, position=pos1, strandedness=True)
        nn1.id = nn1.calculate_id()
        self.assertEqual(new_nodes[1], nn1)
        pos2 = (1008, 1017)
        nn2 = Node(None, block_id=new_bid, path_id=2, position=pos2, strandedness=False)
        nn2.id = nn2.calculate_id()
        self.assertEqual(new_nodes[2], nn2)
        pos3 = (96, 4)
        nn3 = Node(None, block_id=new_bid, path_id=3, position=pos3, strandedness=False)
        nn3.id = nn3.calculate_id()
        self.assertEqual(new_nodes[3], nn3)
        self.assertEqual(new_nodes, {1: nn1, 2: nn2, 3: nn3})

        # block edits
        n1ed = Edit(
            subs=[Substitution(3, "T")],
            dels=[Deletion(8, 2)],
            ins=[Insertion(0, "A")],
        )
        self.assertEqual(new_b.alignment[nn1.id], n1ed)
        n2ed = Edit(
            subs=[Substitution(9, "G")],
            dels=[Deletion(3, 2)],
            ins=[Insertion(7, "T")],
        )
        self.assertEqual(new_b.alignment[nn2.id], n2ed)
        n3ed = Edit(
            subs=[],
            dels=[Deletion(0, 2)],
            ins=[],
        )
        self.assertEqual(new_b.alignment[nn3.id], n3ed)

    def test_block_slice_rev_append(self):
        b, G = self.generate_example()
        new_bid = 42
        i = Interval(
            start=10,
            end=20,
            aligned=True,
            new_block_id=new_bid,
            orientation=False,
            is_anchor=False,
        )

        new_b, new_nodes = block_slice(b, i, G)

        # block consensus
        self.assertEqual(new_b.consensus, "TATATTTATC")

        # new nodes and genome position
        pos1 = (111, 120)
        nn1 = Node(None, block_id=new_bid, path_id=1, position=pos1, strandedness=False)
        nn1.id = nn1.calculate_id()
        self.assertEqual(new_nodes[1], nn1)
        pos2 = (1008, 1017)
        nn2 = Node(None, block_id=new_bid, path_id=2, position=pos2, strandedness=True)
        nn2.id = nn2.calculate_id()
        self.assertEqual(new_nodes[2], nn2)
        pos3 = (96, 4)
        nn3 = Node(None, block_id=new_bid, path_id=3, position=pos3, strandedness=True)
        nn3.id = nn3.calculate_id()
        self.assertEqual(new_nodes[3], nn3)
        self.assertEqual(new_nodes, {1: nn1, 2: nn2, 3: nn3})

        # block edits
        n1ed = Edit(
            subs=[Substitution(3, "T")],
            dels=[Deletion(8, 2)],
            ins=[Insertion(0, "A")],
        )
        self.assertEqual(new_b.alignment[nn1.id], n1ed)
        n2ed = Edit(
            subs=[Substitution(9, "G")],
            dels=[Deletion(3, 2)],
            ins=[Insertion(7, "T")],
        )
        self.assertEqual(new_b.alignment[nn2.id], n2ed)
        n3ed = Edit(
            subs=[],
            dels=[Deletion(0, 2)],
            ins=[],
        )
        self.assertEqual(new_b.alignment[nn3.id], n3ed)

    def test_sequence_reconstruction(self):
        b, G = self.generate_example()
        s1 = apply_edits_to_ref(b.alignment[1], b.consensus)
        s2 = apply_edits_to_ref(b.alignment[2], b.consensus)
        s3 = apply_edits_to_ref(b.alignment[3], b.consensus)
        i1 = Interval(
            start=0,
            end=10,
            aligned=True,
            new_block_id=2,
            orientation=True,
            is_anchor=True,
        )
        i2 = Interval(
            start=10,
            end=20,
            aligned=True,
            new_block_id=2,
            orientation=True,
            is_anchor=True,
        )
        i3 = Interval(
            start=20,
            end=27,
            aligned=True,
            new_block_id=2,
            orientation=True,
            is_anchor=True,
        )
        nb1, nn1 = block_slice(b, i1, G)
        nb2, nn2 = block_slice(b, i2, G)
        nb3, nn3 = block_slice(b, i3, G)

        for i_old, s_old in [(1, s1), (2, s2), (3, s3)]:
            s_new_1 = apply_edits_to_ref(nb1.alignment[nn1[i_old].id], nb1.consensus)
            s_new_2 = apply_edits_to_ref(nb2.alignment[nn2[i_old].id], nb2.consensus)
            s_new_3 = apply_edits_to_ref(nb3.alignment[nn3[i_old].id], nb3.consensus)

            self.assertEqual(s_old, s_new_1 + s_new_2 + s_new_3)


if __name__ == "__main__":
    unittest.main()
