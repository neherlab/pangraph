from functions import *
import unittest


class TestGraphUpdate(unittest.TestCase):
    def test_graph_update(self):
        # graph
        # p1 -> [b1+,b2+,b3+]
        # p2 -> [b2+,b3+]
        # p3 -> [b1+,b2-,b3+]
        # update
        # b2+ -> [b4+, b5-]

        # nodes
        n1_b1_p1 = Node(id=1, block_id=1, path_id=1, strandedness=True, position=())
        n2_b1_p3 = Node(id=2, block_id=1, path_id=3, strandedness=True, position=())
        n3_b2_p1 = Node(id=3, block_id=2, path_id=1, strandedness=True, position=())
        n4_b2_p2 = Node(id=4, block_id=2, path_id=2, strandedness=True, position=())
        n5_b2_p3 = Node(id=5, block_id=2, path_id=3, strandedness=False, position=())
        n6_b3_p1 = Node(id=6, block_id=3, path_id=1, strandedness=True, position=())
        n7_b3_p2 = Node(id=7, block_id=3, path_id=2, strandedness=True, position=())
        n8_b3_p3 = Node(id=8, block_id=3, path_id=3, strandedness=True, position=())
        nodes = {
            1: n1_b1_p1,
            2: n2_b1_p3,
            3: n3_b2_p1,
            4: n4_b2_p2,
            5: n5_b2_p3,
            6: n6_b3_p1,
            7: n7_b3_p2,
            8: n8_b3_p3,
        }

        # blocks
        b1 = Block(id=1, consensus="1", alignment={1: None, 2: None})
        b2 = Block(id=2, consensus="2", alignment={3: None, 4: None, 5: None})
        b3 = Block(id=3, consensus="3", alignment={6: None, 7: None, 8: None})
        blocks = {1: b1, 2: b2, 3: b3}

        # paths
        p1 = Path(id=1, nodes=[1, 3, 6], L=None)
        p2 = Path(id=2, nodes=[4, 7], L=None)
        p3 = Path(id=3, nodes=[2, 5, 8], L=None)
        paths = {1: p1, 2: p2, 3: p3}

        # graph
        G = Pangraph(paths=paths, blocks=blocks, nodes=nodes)

        # new nodes
        n9_b4_p1 = Node(id=9, block_id=4, path_id=1, strandedness=True, position=())
        n10_b5_p1 = Node(id=10, block_id=5, path_id=1, strandedness=False, position=())
        n11_b4_p2 = Node(id=11, block_id=4, path_id=2, strandedness=True, position=())
        n12_b5_p2 = Node(id=12, block_id=5, path_id=2, strandedness=False, position=())
        n13_b4_p3 = Node(id=13, block_id=4, path_id=3, strandedness=False, position=())
        n14_b5_p3 = Node(id=14, block_id=5, path_id=3, strandedness=True, position=())
        # new blocks
        b4 = Block(id=4, consensus="4", alignment={9: None, 11: None, 13: None})
        b5 = Block(id=5, consensus="5", alignment={10: None, 12: None, 14: None})
        # update
        u = GraphUpdate(
            b_old_id=2,
            b_new=[b4, b5],
            n_new={
                3: [n9_b4_p1, n10_b5_p1],
                4: [n11_b4_p2, n12_b5_p2],
                5: [n14_b5_p3, n13_b4_p3],
            },
        )

        # apply update
        G.update(u)

        # check blocks
        self.assertEqual(G.blocks, {1: b1, 3: b3, 4: b4, 5: b5})
        # check paths
        p1_new = Path(id=1, nodes=[1, 9, 10, 6], L=None)
        p2_new = Path(id=2, nodes=[11, 12, 7], L=None)
        p3_new = Path(id=3, nodes=[2, 14, 13, 8], L=None)
        self.assertEqual(G.paths, {1: p1_new, 2: p2_new, 3: p3_new})
        # check nodes
        new_nodes = {
            1: n1_b1_p1,
            2: n2_b1_p3,
            6: n6_b3_p1,
            7: n7_b3_p2,
            8: n8_b3_p3,
            9: n9_b4_p1,
            10: n10_b5_p1,
            11: n11_b4_p2,
            12: n12_b5_p2,
            13: n13_b4_p3,
            14: n14_b5_p3,
        }
        self.assertEqual(G.nodes, new_nodes)


class TestTargetBlocks(unittest.TestCase):
    def test_target_blocks(self):
        new_hit = lambda name: Hit(name=name, length=None, start=None, stop=None)
        new_aln = lambda qry, reff: Alignment(
            qry=qry,
            reff=reff,
            orientation=None,
        )

        h1 = new_hit(1)
        h2 = new_hit(2)
        h3 = new_hit(3)
        h4 = new_hit(4)
        h5 = new_hit(1)
        h6 = new_hit(2)
        h7 = new_hit(3)
        h8 = new_hit(4)

        a1 = new_aln(h1, h2)  # a1 : 1 -- 2
        a2 = new_aln(h3, h4)  # a2 : 3 -- 4
        a3 = new_aln(h5, h8)  # a3 : 1 -- 4
        a4 = new_aln(h6, h7)  # a4 : 2 -- 3

        TB = target_blocks([a1, a2, a3, a4])

        self.assertEqual(TB, {1: [a1, a3], 2: [a1, a4], 3: [a2, a4], 4: [a2, a3]})


class TestExtractHits(unittest.TestCase):
    def test_extract_hits(self):
        bid = 1
        new_hit = lambda name, start: Hit(
            name=name, length=None, start=start, stop=None
        )
        create_hit = lambda new_bid, deep, strand, hit: {
            "new_block_id": new_bid,
            "deep": deep,
            "orientation": strand,
            "hit": hit,
        }
        new_aln = lambda id, reff, qry, strand, deep: Alignment(
            reff=reff,
            qry=qry,
            new_block_id=id,
            orientation=strand,
            deep_block=deep,
        )

        h1_a = new_hit(1, 10)
        h1_b = new_hit(1, 20)
        h1_c = new_hit(1, 30)
        h1_d = new_hit(1, 40)
        h2_e = new_hit(2, 50)
        h2_f = new_hit(2, 60)
        h2_g = new_hit(2, 70)
        h2_h = new_hit(2, 80)

        a1 = new_aln(3, h1_a, h1_b, True, "reff")  # 1 -- 1
        a2 = new_aln(4, h1_c, h2_e, True, "qry")  # 1 -- 2
        a3 = new_aln(5, h2_f, h1_d, False, "reff")  # 2 -- 1
        a4 = new_aln(6, h2_g, h2_h, False, "qry")  # 2 -- 2

        hits = extract_hits(bid, [a1, a2, a3, a4])

        self.assertEqual(
            hits,
            [
                create_hit(3, True, True, h1_a),
                create_hit(3, False, True, h1_b),
                create_hit(4, False, True, h1_c),
                create_hit(5, False, False, h1_d),
            ],
        )


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
