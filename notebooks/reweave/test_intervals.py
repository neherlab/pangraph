from intervals_utils import *
import unittest


class TestIntervals(unittest.TestCase):
    def example(self):
        block_L = 1000
        bid = 0
        create_hit = lambda new_bid, deep, strand, start, stop: {
            "new_block_id": new_bid,
            "deep": deep,
            "orientation": strand,
            "hit": Hit(name=bid, length=None, start=start, stop=stop),
        }
        hits = [
            create_hit(1, True, True, 10, 100),
            create_hit(2, False, True, 200, 300),
            create_hit(3, True, True, 310, 500),
            create_hit(4, False, True, 600, 900),
        ]
        return hits, block_L, bid

    def test_create_intervals(self):
        hits, block_L, bid = self.example()
        I = create_intervals(hits, block_L)
        self.assertEqual(
            I,
            [
                Interval(
                    start=0, end=10, aligned=False, new_block_id=hash((bid, 0, 10))
                ),
                Interval(
                    start=10,
                    end=100,
                    aligned=True,
                    new_block_id=1,
                    deep=True,
                    orientation=True,
                ),
                Interval(
                    start=100,
                    end=200,
                    aligned=False,
                    new_block_id=hash((bid, 100, 200)),
                ),
                Interval(
                    start=200,
                    end=300,
                    aligned=True,
                    new_block_id=2,
                    deep=False,
                    orientation=True,
                ),
                Interval(
                    start=300,
                    end=310,
                    aligned=False,
                    new_block_id=hash((bid, 300, 310)),
                ),
                Interval(
                    start=310,
                    end=500,
                    aligned=True,
                    new_block_id=3,
                    deep=True,
                    orientation=True,
                ),
                Interval(
                    start=500,
                    end=600,
                    aligned=False,
                    new_block_id=hash((bid, 500, 600)),
                ),
                Interval(
                    start=600,
                    end=900,
                    aligned=True,
                    new_block_id=4,
                    deep=False,
                    orientation=True,
                ),
                Interval(
                    start=900,
                    end=1000,
                    aligned=False,
                    new_block_id=hash((bid, 900, 1000)),
                ),
            ],
        )

    def test_refine_intervals(self):
        hits, block_L, bid = self.example()
        thr_len = 50
        I = extract_intervals(hits, block_L, thr_len=thr_len)
        self.assertEqual(
            I,
            [
                Interval(
                    start=0,
                    end=100,
                    aligned=True,
                    new_block_id=1,
                    deep=True,
                    orientation=True,
                ),
                Interval(
                    start=100,
                    end=200,
                    aligned=False,
                    new_block_id=hash((bid, 100, 200)),
                ),
                Interval(
                    start=200,
                    end=300,
                    aligned=True,
                    new_block_id=2,
                    deep=False,
                    orientation=True,
                ),
                Interval(
                    start=300,
                    end=500,
                    aligned=True,
                    new_block_id=3,
                    deep=True,
                    orientation=True,
                ),
                Interval(
                    start=500,
                    end=600,
                    aligned=False,
                    new_block_id=hash((bid, 500, 600)),
                ),
                Interval(
                    start=600,
                    end=900,
                    aligned=True,
                    new_block_id=4,
                    deep=False,
                    orientation=True,
                ),
                Interval(
                    start=900,
                    end=1000,
                    aligned=False,
                    new_block_id=hash((bid, 900, 1000)),
                ),
            ],
        )


if __name__ == "__main__":
    unittest.main()
