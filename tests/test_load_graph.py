import pypangraph as pp


def test_load_graph():
    fname = "tests/data/plasmid_graph.json"
    pan = pp.Pangraph.load_json(fname)
