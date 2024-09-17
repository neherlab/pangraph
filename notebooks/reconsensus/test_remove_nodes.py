import pytest

from utils import Edit, Insertion, Deletion, Substitution
from copy import deepcopy

from reconsensus import (
    Block,
    Node,
    Path,
    Pangraph,
    find_empty_nodes,
    remove_emtpy_nodes,
)


@pytest.fixture
def input_graph():
    nodes = {
        1: Node(1, 1, 1, (0, 10), True),
        2: Node(2, 1, 2, (0, 10), True),
        3: Node(3, 1, 3, (0, 0), False),
        4: Node(4, 2, 1, (10, 20), True),
        5: Node(5, 2, 3, (0, 10), True),
    }
    paths = {
        1: Path(1, [1, 4], 20, False),
        2: Path(2, [2], 10, False),
        3: Path(3, [3, 5], 10, False),
    }
    blocks = {
        1: Block(
            1,
            "AAAAAAAAAA",
            alignment={
                1: Edit(ins=[Insertion(0, "GGG")], dels=[Deletion(1, 3)], subs=[]),
                2: Edit(ins=[], dels=[], subs=[Substitution(5, "G")]),
                3: Edit(ins=[], dels=[Deletion(0, 10)], subs=[]),
            },
        ),
        2: Block(
            2,
            "CCCCCCCCCC",
            alignment={
                4: Edit(ins=[], dels=[], subs=[]),
                5: Edit(ins=[], dels=[], subs=[]),
            },
        ),
    }
    graph = Pangraph(paths, blocks, nodes)
    return graph


@pytest.fixture
def graph_remove_node(input_graph):
    graph = deepcopy(input_graph)
    remove_emtpy_nodes(graph, block_ids=[1, 2])
    return graph


def test_find_empty_nodes(input_graph):
    empty_node_ids = find_empty_nodes(input_graph, [1, 2])
    assert empty_node_ids == [3]


def test_removed_node(graph_remove_node):
    assert set(graph_remove_node.nodes.keys()) == {1, 2, 4, 5}


def test_cleaned_paths(graph_remove_node):
    assert graph_remove_node.paths[3].nodes == [5]


def test_cleaned_alignments(graph_remove_node):
    assert 3 not in graph_remove_node.blocks[1].alignment


def test_same_nodes(graph_remove_node, input_graph):
    del input_graph.nodes[3]
    assert input_graph.nodes == graph_remove_node.nodes


def test_same_paths(graph_remove_node, input_graph):
    input_graph.paths[3].nodes.remove(3)
    assert input_graph.paths == graph_remove_node.paths


def test_same_blocks(graph_remove_node, input_graph):
    del input_graph.blocks[1].alignment[3]
    assert input_graph.blocks == graph_remove_node.blocks
