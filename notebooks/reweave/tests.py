from utils import *
from test_slice import *
from test_intervals import *
import unittest
import numpy as np
from copy import deepcopy


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


class TestExtractHits(unittest.TestCase):
    def test_extract_hits(self):
        bid = 1
        new_hit = lambda name, start: Hit(
            name=name, length=None, start=start, stop=None
        )
        create_hit = lambda new_bid, is_anchor, strand, hit: {
            "new_block_id": new_bid,
            "is_anchor": is_anchor,
            "orientation": strand,
            "hit": hit,
        }
        new_aln = lambda id, reff, qry, strand, anc: Alignment(
            reff=reff,
            qry=qry,
            new_block_id=id,
            orientation=strand,
            anchor_block=anc,
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


class TestGroupPromises(unittest.TestCase):
    def test_group_promises(self):
        b1_anchor = Block(id=1, consensus="A", alignment={1: None, 2: None, 3: None})
        b1_append = Block(id=1, consensus="B", alignment={4: None, 5: None})
        b2_anchor = Block(id=2, consensus="C", alignment={6: None, 7: None, 8: None})
        b2_append = Block(id=2, consensus="D", alignment={7: None, 8: None})
        b3_anchor = Block(id=3, consensus="E", alignment={11: None, 12: None})
        b3_append = Block(id=3, consensus="F", alignment={13: None})

        t = lambda b, a, o: ToMerge(block=b, orientation=o, is_anchor=a)

        H = [
            t(b1_anchor, True, True),
            t(b1_append, False, True),
            t(b3_anchor, True, False),
            t(b2_append, False, True),
            t(b2_anchor, True, True),
            t(b3_append, False, False),
        ]
        promises = group_promises(H)
        self.assertEqual(
            promises,
            [
                MergePromise(
                    anchor_block=b1_anchor, append_block=b1_append, orientation=True
                ),
                MergePromise(
                    anchor_block=b2_anchor, append_block=b2_append, orientation=True
                ),
                MergePromise(
                    anchor_block=b3_anchor, append_block=b3_append, orientation=False
                ),
            ],
        )


class TestMisc(unittest.TestCase):

    def test_assign_anchor_block(self):
        b1 = Block(id=1, consensus="A", alignment={1: None, 2: None, 3: None})
        b2 = Block(id=2, consensus="B", alignment={4: None, 5: None})
        b3 = Block(id=3, consensus="C", alignment={6: None})
        b4 = Block(id=4, consensus="D", alignment={7: None, 8: None, 9: None, 10: None})
        G = Pangraph(paths={}, blocks={1: b1, 2: b2, 3: b3, 4: b4}, nodes={})

        h = lambda name: Hit(name=name, length=None, start=None, stop=None)
        a = lambda q, r: Alignment(qry=h(q), reff=h(r), orientation=None)
        mergers = [
            a(1, 2),
            a(3, 4),
            a(4, 1),
        ]
        assign_anchor_block(mergers, G)

        self.assertEqual(mergers[0].anchor_block, "qry")
        self.assertEqual(mergers[1].anchor_block, "reff")
        self.assertEqual(mergers[2].anchor_block, "qry")

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


class TestSplitBlock(unittest.TestCase):

    def generate_example(self):
        bid = 1
        nid1 = 1000
        nid2 = 2000
        nid3 = 3000

        # nodes
        n1 = Node(
            id=nid1, block_id=bid, path_id=100, strandedness=True, position=(100, 230)
        )
        n2 = Node(
            id=2000,
            block_id=bid,
            path_id=200,
            strandedness=False,
            position=(1000, 1130),
        )
        n3 = Node(
            id=nid2, block_id=bid, path_id=300, strandedness=False, position=(180, 110)
        )

        # block
        np.random.seed(1)
        gseq = lambda l: "".join([np.random.choice(list("ACGT")) for _ in range(l)])
        ed1 = Edit(ins=[], dels=[], subs=[])
        ed2 = Edit(ins=[], dels=[], subs=[])
        ed3 = Edit(ins=[], dels=[], subs=[])
        b1 = Block(
            id=bid, consensus=gseq(130), alignment={nid1: ed1, nid2: ed2, nid3: ed3}
        )

        # paths
        p1 = Path(id=100, nodes=[nid1], L=2000)
        p2 = Path(id=200, nodes=[nid2], L=2000)
        p3 = Path(id=300, nodes=[nid3], L=200)
        G = Pangraph(
            paths={100: p1, 200: p2, 300: p3},
            blocks={bid: b1},
            nodes={nid1: n1, nid2: n2, nid3: n3},
        )

        # alignments
        h = lambda name, start, stop: Hit(
            name=name, length=None, start=start, stop=stop
        )
        a = lambda q, r, strand, nid, anc: Alignment(
            qry=q, reff=r, orientation=strand, new_block_id=nid, anchor_block=anc
        )
        M = [
            a(h(1, 10, 50), h(2, 100, 150), True, 42, "qry"),
            a(h(3, 1000, 1050), h(1, 80, 130), False, 43, "qry"),
        ]

        return G, M, bid

    def test_split_block(self):

        thr_len = 20
        G, M, bid = self.generate_example()
        u, H = split_block(bid, M, G, thr_len)

        # check graph update
        self.assertEqual(u.b_old_id, bid)

        # check node updates
        nodekeys_1 = []
        nodekeys_2 = []
        nodekeys_3 = []
        n1, n2, n3 = u.n_new[1000]
        nodekeys_1.append(n1.id)
        nodekeys_2.append(n2.id)
        nodekeys_3.append(n3.id)
        self.assertEqual(n1.block_id, 42)
        self.assertEqual(n1.strandedness, True)
        self.assertEqual(n1.position, (100, 150))
        self.assertEqual(n2.strandedness, True)
        self.assertEqual(n2.position, (150, 180))
        self.assertEqual(n3.block_id, 43)
        self.assertEqual(n3.strandedness, False)
        self.assertEqual(n3.position, (180, 230))

        n1, n2, n3 = u.n_new[2000]
        nodekeys_1.append(n3.id)
        nodekeys_2.append(n2.id)
        nodekeys_3.append(n1.id)
        self.assertEqual(n1.block_id, 43)
        self.assertEqual(n1.strandedness, True)
        self.assertEqual(n1.position, (1000, 1050))
        self.assertEqual(n2.strandedness, False)
        self.assertEqual(n2.position, (1050, 1080))
        self.assertEqual(n3.block_id, 42)
        self.assertEqual(n3.strandedness, False)
        self.assertEqual(n3.position, (1080, 1130))

        n1, n2, n3 = u.n_new[3000]
        nodekeys_1.append(n3.id)
        nodekeys_2.append(n2.id)
        nodekeys_3.append(n1.id)
        self.assertEqual(n1.block_id, 43)
        self.assertEqual(n1.strandedness, True)
        self.assertEqual(n1.position, (180, 30))
        self.assertEqual(n2.strandedness, False)
        self.assertEqual(n2.position, (30, 60))
        self.assertEqual(n3.block_id, 42)
        self.assertEqual(n3.strandedness, False)
        self.assertEqual(n3.position, (60, 110))

        # check block updates
        self.assertEqual(len(u.b_new), 1)
        b = u.b_new[0]
        self.assertEqual(b.consensus, G.blocks[bid].consensus[50:80])
        self.assertSetEqual(set(b.alignment.keys()), set(nodekeys_2))

        # check ToMerge objects
        self.assertEqual(len(H), 2)
        H1, H2 = H
        self.assertEqual(H1.block.id, 42)
        self.assertEqual(H1.is_anchor, True)
        self.assertEqual(H1.orientation, True)
        self.assertEqual(H2.block.id, 43)
        self.assertEqual(H2.is_anchor, False)
        self.assertEqual(H2.orientation, False)
        b1, b3 = H1.block, H2.block
        self.assertEqual(b1.consensus, G.blocks[bid].consensus[0:50])
        self.assertEqual(b3.consensus, G.blocks[bid].consensus[80:130])
        self.assertEqual(set(b1.alignment.keys()), set(nodekeys_1))
        self.assertEqual(set(b3.alignment.keys()), set(nodekeys_3))


class TestReweave(unittest.TestCase):

    def generate_example(self):
        np.random.seed(0)

        # block lengths
        Bl = {
            10: 200,
            20: 400,
            30: 100,
            40: 500,
            50: 250,
        }
        Bseq = {k: "".join(np.random.choice(list("ACGT"), l)) for k, l in Bl.items()}

        # nodes
        n = lambda idx, bid, path, pos, s: Node(idx, bid, path, pos, s)
        nodes = {
            1: n(1, 10, 100, (700, 885), True),
            2: n(2, 30, 100, (885, 988), True),
            3: n(3, 30, 200, (100, 180), False),
            4: n(4, 20, 200, (180, 555), False),
            5: n(5, 10, 200, (555, 735), False),
            6: n(6, 40, 300, (600, 100), True),
            7: n(7, 50, 300, (100, 325), True),
            8: n(8, 50, 300, (325, 580), False),
        }

        # paths
        paths = {
            100: Path(100, [1, 2], 1000),
            200: Path(200, [3, 4, 5], 1000),
            300: Path(300, [6, 7, 8], 1000),
        }

        # edits
        i = lambda p, l, a: Insertion(p, a * l)
        d = lambda p, l: Deletion(p, l)
        s = lambda p, a: Substitution(p, a)
        ed = {
            1: Edit(ins=[i(150, 10, "T")], dels=[d(50, 25)], subs=[s(125, "G")]),
            2: Edit(ins=[i(50, 3, "G")], dels=[], subs=[]),
            3: Edit(ins=[i(25, 5, "G")], dels=[d(50, 25)], subs=[]),
            4: Edit(
                ins=[i(250, 5, "A"), i(300, 5, "A")],
                dels=[d(100, 25), d(350, 10)],
                subs=[s(50, "G"), s(225, "T")],
            ),
            5: Edit(ins=[i(200, 5, "A")], dels=[d(100, 25)], subs=[s(25, "T")]),
            6: Edit(ins=[i(200, 10, "T")], dels=[d(350, 10)], subs=[s(100, "T")]),
            7: Edit(ins=[], dels=[d(100, 25)], subs=[s(50, "G")]),
            8: Edit(ins=[i(150, 5, "T")], dels=[], subs=[]),
        }

        # blocks
        b = lambda bid, nids: Block(bid, Bseq[bid], {nid: ed[nid] for nid in nids})
        blocks = {
            10: b(10, [1, 5]),
            20: b(20, [4]),
            30: b(30, [2, 3]),
            40: b(40, [6]),
            50: b(50, [7, 8]),
        }

        # graph
        G = Pangraph(paths, blocks, nodes)

        # alignments
        M = [
            Alignment(
                qry=Hit(10, 200, 10, 200),  # anchor
                reff=Hit(40, 500, 10, 200),
                orientation=True,
            ),
            Alignment(
                qry=Hit(20, 400, 0, 200),
                reff=Hit(40, 500, 300, 500),  # anchor
                orientation=False,
            ),
            Alignment(
                qry=Hit(20, 400, 300, 400),
                reff=Hit(50, 250, 0, 100),  # anchor
                orientation=True,
            ),
            Alignment(
                qry=Hit(30, 100, 0, 100),
                reff=Hit(50, 250, 150, 250),  # anchor
                orientation=True,
            ),
        ]

        return G, M

    def test_reweave(self):
        i = lambda p, l, a: Insertion(p, a * l)
        d = lambda p, l: Deletion(p, l)
        s = lambda p, a: Substitution(p, a)

        G, M = self.generate_example()
        O = deepcopy(G)
        thr_len = 90
        G, P = reweave(M, G, thr_len=90)

        # new paths
        p1 = G.paths[100]
        p2 = G.paths[200]
        p3 = G.paths[300]

        # new nodes
        self.assertEqual(len(p1.nodes), 2)
        nid_100_1, nid_100_2 = p1.nodes
        self.assertEqual(len(p2.nodes), 5)
        nid_200_1, nid_200_2, nid_200_3, nid_200_4, nid_200_5 = p2.nodes
        self.assertEqual(len(p3.nodes), 7)
        nid_300_1, nid_300_2, nid_300_3, nid_300_4, nid_300_5, nid_300_6, nid_300_7 = (
            p3.nodes
        )

        # check node positions
        self.assertEqual(G.nodes[nid_100_1].position, O.nodes[1].position)
        self.assertEqual(G.nodes[nid_100_2].position, O.nodes[2].position)
        self.assertEqual(G.nodes[nid_200_1].position, O.nodes[3].position)
        self.assertEqual(G.nodes[nid_200_2].position, (180, 275))
        self.assertEqual(G.nodes[nid_200_3].position, (275, 380))
        self.assertEqual(G.nodes[nid_200_4].position, (380, 555))
        self.assertEqual(G.nodes[nid_200_5].position, O.nodes[5].position)
        self.assertEqual(G.nodes[nid_300_1].position, (600, 800))
        self.assertEqual(G.nodes[nid_300_2].position, (800, 910))
        self.assertEqual(G.nodes[nid_300_3].position, (910, 100))
        self.assertEqual(G.nodes[nid_300_4].position, (100, 225))
        self.assertEqual(G.nodes[nid_300_5].position, (225, 325))
        self.assertEqual(G.nodes[nid_300_6].position, (325, 430))
        self.assertEqual(G.nodes[nid_300_7].position, (430, 580))

        # check node orientation
        self.assertEqual(G.nodes[nid_100_1].strandedness, True)
        self.assertEqual(G.nodes[nid_100_2].strandedness, True)
        self.assertEqual(G.nodes[nid_200_1].strandedness, False)
        self.assertEqual(G.nodes[nid_200_2].strandedness, False)
        self.assertEqual(G.nodes[nid_200_3].strandedness, False)
        self.assertEqual(G.nodes[nid_200_4].strandedness, True)
        self.assertEqual(G.nodes[nid_200_5].strandedness, False)
        self.assertEqual(G.nodes[nid_300_1].strandedness, True)
        self.assertEqual(G.nodes[nid_300_2].strandedness, True)
        self.assertEqual(G.nodes[nid_300_3].strandedness, True)
        self.assertEqual(G.nodes[nid_300_4].strandedness, True)
        self.assertEqual(G.nodes[nid_300_5].strandedness, True)
        self.assertEqual(G.nodes[nid_300_6].strandedness, False)
        self.assertEqual(G.nodes[nid_300_7].strandedness, False)

        # block identity
        bid10_1 = G.nodes[nid_100_1].block_id
        self.assertEqual(G.nodes[nid_200_5].block_id, bid10_1)
        self.assertNotIn(bid10_1, G.blocks)  # missing block key, it is in merge promise
        self.assertIn(bid10_1, [p.anchor_block.id for p in P])

        bid20_2 = G.nodes[nid_200_3].block_id
        self.assertIn(bid20_2, G.blocks)  # block is already in the graph
        self.assertNotIn(
            bid20_2, [p.anchor_block.id for p in P] + [p.append_block.id for p in P]
        )
        ed20_2 = G.blocks[bid20_2].alignment[nid_200_3]
        self.assertEqual(ed20_2, Edit(ins=[i(50, 5, "A")], dels=[], subs=[s(25, "T")]))

        bid20_1 = G.nodes[nid_200_1].block_id
        for n in [nid_100_2, nid_300_5, nid_300_6]:
            self.assertEqual(G.nodes[n].block_id, bid20_1)

        # merge promises
        self.assertEqual(len(P), 4)
        P_dict = {p.anchor_block.id: p for p in P}

        p1 = P_dict[bid10_1]
        self.assertEqual(p1.orientation, True)
        self.assertEqual(p1.anchor_block.consensus, O.blocks[10].consensus)
        self.assertEqual(p1.append_block.consensus, O.blocks[40].consensus[0:200])
        self.assertEqual(p1.append_block.id, G.nodes[nid_300_1].block_id)

        bid40_3 = G.nodes[nid_300_3].block_id
        p2 = P_dict[bid40_3]
        self.assertEqual(p2.orientation, False)
        self.assertEqual(p2.anchor_block.consensus, O.blocks[40].consensus[300:500])
        self.assertEqual(p2.append_block.id, G.nodes[nid_200_4].block_id)
        self.assertEqual(p2.append_block.consensus, O.blocks[20].consensus[0:200])


if __name__ == "__main__":
    unittest.main()
