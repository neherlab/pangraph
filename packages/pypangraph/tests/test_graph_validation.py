"""Schema-validation behaviour of ``Pangraph.from_json``.

Before constructing the object, the loader validates every graph against the JSON
schema that is generated from the Rust types. These tests lock in which graphs are
accepted and which are rejected, so the contract holds no matter which validation
engine backs the loader.

Oracle: ``pypangraph/pangraph_schema.py`` (the schema itself). Each rejection case
violates one named schema constraint.
"""

import json

import pytest

import pypangraph as pp


def _write(tmp_path, graph):
    fname = tmp_path / "graph.json"
    fname.write_text(json.dumps(graph))
    return fname


def test_graph_validation_accepts_valid(tmp_path, junction_pangraph_json):
    fname = _write(tmp_path, junction_pangraph_json)
    pan = pp.Pangraph.from_json(fname)
    assert len(pan.strains()) == 3


def _drop_blocks(graph):
    # Pangraph.blocks is a required property.
    del graph["blocks"]


def _nodes_as_list(graph):
    # Pangraph.nodes must be an object (map), not an array.
    graph["nodes"] = [1, 2, 3]


def _negative_tot_len(graph):
    # PangraphPath.tot_len is a uint (minimum 0).
    next(iter(graph["paths"].values()))["tot_len"] = -5


def _bad_strand(graph):
    # Strand is the enum {"+", "-"}.
    next(iter(graph["nodes"].values()))["strand"] = "x"


def _missing_node_block_id(graph):
    # PangraphNode.block_id is required.
    del next(iter(graph["nodes"].values()))["block_id"]


def _consensus_not_string(graph):
    # PangraphBlock.consensus must be a string.
    next(iter(graph["blocks"].values()))["consensus"] = 123


@pytest.mark.parametrize(
    "mutate",
    [
        pytest.param(_drop_blocks, id="missing_blocks"),
        pytest.param(_nodes_as_list, id="nodes_as_list"),
        pytest.param(_negative_tot_len, id="negative_tot_len"),
        pytest.param(_bad_strand, id="bad_strand_enum"),
        pytest.param(_missing_node_block_id, id="missing_node_block_id"),
        pytest.param(_consensus_not_string, id="consensus_not_string"),
    ],
)
def test_graph_validation_rejects_malformed(tmp_path, junction_pangraph_json, mutate):
    mutate(junction_pangraph_json)
    fname = _write(tmp_path, junction_pangraph_json)
    with pytest.raises(pp.PangraphLoadError, match="invalid pangraph JSON"):
        pp.Pangraph.from_json(fname)
