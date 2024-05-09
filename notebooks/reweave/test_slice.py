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


if __name__ == "__main__":
    unittest.main()
