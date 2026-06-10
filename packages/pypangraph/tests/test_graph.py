import pytest
import pypangraph as pp


def test_load_graph():
    fname = "tests/data/plasmids.json"
    pp.Pangraph.from_json(fname)


@pytest.fixture
def graph():
    fname = "tests/data/plasmids.json"
    pan = pp.Pangraph.from_json(fname)
    return pan


def test_load_graph_gz():
    fname = "tests/data/plasmids.json.gz"
    pan = pp.Pangraph.from_json(fname)
    assert len(pan.strains()) == 15


def test_load_graph_invalid_extension(tmp_path):
    fname = tmp_path / "plasmids.txt"
    fname.write_text("{}")

    with pytest.raises(pp.PangraphLoadError, match=r"\.json or \.json\.gz"):
        pp.Pangraph.from_json(fname)


def test_load_graph_invalid_json(tmp_path):
    fname = tmp_path / "broken.json"
    fname.write_text("{this is not valid json}")

    with pytest.raises(pp.PangraphLoadError, match="failed to load pangraph"):
        pp.Pangraph.from_json(fname)


def test_load_graph_invalid_schema(tmp_path):
    fname = tmp_path / "invalid.json"
    fname.write_text("{}")

    with pytest.raises(pp.PangraphLoadError, match="invalid pangraph JSON"):
        pp.Pangraph.from_json(fname)


def test_paths(graph):
    path = graph.paths["RCS48_p1"]
    assert len(path) == 60
    assert path.nuc_len == 80596
    assert graph.paths.idx_to_name[0] == "RCS48_p1"


def test_get_strains(graph):
    strains = graph.strains()
    assert len(strains) == 15


def test_blockstats_df(graph):
    df = graph.to_blockstats_df()
    assert df.shape[0] == 137
    assert df["core"].sum() == 27
    assert df["duplicated"].sum() == 10


def test_block_ids_are_strings(graph):
    """Block ids are stored as strings internally, and the blockstats index is str.

    Block ids are u64 hashes that overflow int64; storing them as strings keeps them
    out of any float64 coercion (which would silently corrupt values above 2**53).
    """
    bdf = graph.to_blockstats_df()
    assert bdf.index.dtype == object
    assert all(isinstance(bid, str) for bid in bdf.index)
    assert all(isinstance(bid, str) for bid in graph.blocks.keys())


def test_blocks_accessor_accepts_int_or_str(graph):
    """pan.blocks[bid] resolves the same block whether bid is passed as str or int."""
    bid = graph.blocks.keys()[0]  # str
    assert isinstance(bid, str)
    assert graph.blocks[int(bid)] is graph.blocks[bid]
    assert int(bid) in graph.blocks and bid in graph.blocks


def test_blockcount_df(graph):
    df = graph.to_blockcount_df()
    assert df.shape[0] == 137
    assert df.shape[1] == 15
    assert df.sum().sum() == 1042


def test_nodes_to_blocks(graph):
    path = graph.paths["RCS49_p1"]
    nodes = path.nodes
    B, S = graph.nodes.nodes_to_blocks(nodes)
    assert len(B) == len(nodes)
    assert len(S) == len(nodes)

    b, s = graph.nodes.node_to_block(8533989107945450583)
    assert b == "14710008249239879492"  # block ids are strings internally
    assert s


def test_core_genome_alignment(graph):
    core_aln = graph.core_genome_alignment()
    assert len(core_aln) == 15
    assert core_aln.get_alignment_length() == 64989


def test_core_genome_alignment_invalid_guide_strain(graph):
    """An unknown guide strain is rejected with a ValueError (not a bare assert)."""
    with pytest.raises(ValueError, match="Guide strain .* not found"):
        graph.core_genome_alignment(guide_strain="does_not_exist")


def test_pairwise_accessory_genome_comparisons(graph):
    ddf = graph.pairwise_accessory_genome_comparison()
    assert ddf.shape == (225, 2)
    assert ddf.loc["RCS48_p1", "RCS48_p1"]["diff"] == 0
    assert ddf.loc["RCS48_p1", "RCS48_p1"]["shared"] == 79580
    assert ddf.loc["RCS48_p1", "RCS49_p1"]["diff"] == 689
    assert ddf.loc["RCS48_p1", "RCS49_p1"]["shared"] == 79249
