---
sidebar_position: 3
---

# Loading and exploring a graph

We start the tutorial by loading a pangenome graph object and exploring its properties. For this tutorial we will use the `plasmids.json` file that you can find in pypangraph repository under `tests/data/plasmids.json`

We can load the graph object with:

```python
import pypangraph as pp

graph = pp.Pangraph.from_json("plasmids.json")
print(graph)
# pangraph object with 15 paths, 137 blocks and 1042 nodes
```

## The components of a graph

As explained in the [tutorial](../tutorial/t01-building-pangraph.md#what-is-a-pangraph), a pangenome graph is composed of three main components: nodes, blocks and paths.

- **Blocks** encode multiple sequence alignments that group together homologous parts of the input genomes.
- **Paths** are representation of the input genomes as a sequences of blocks. More precisely, as sequence of **nodes**.
- **Nodes** are the elements that connect blocks and paths. Each node connects a single element in a path with a single entry in a block alignment.

These are saved in three sub-properties of the graph object:

```python
graph.blocks
graph.paths
graph.nodes
```

### Paths

This is a graph made from 15 different plasmid sequences. Each sequence corresponds to a different path. We can get the path names with:

```python
print(graph.paths.keys())
# ['RCS48_p1', 'RCS49_p1', 'RCS64_p2', 'RCS80_p1', ... ]
```

Paths are also numbered, and the connection between the name and the number can be retrieved with:

```python
print(graph.paths.idx_to_name)
# {0: 'RCS48_p1', 1: 'RCS49_p1', 2: 'RCS64_p2', 3: 'RCS80_p1', ... }
```

We can recover a specific path with its identifier:

```python
path = graph.paths["RCS48_p1"]
print(path)
# path object | name = RCS48_p1, n. nodes = 60, length = 80596 bp
```

And get the list of node ids that compose the path:

```python
print(path.nodes)
# [11788816242159313242, 6289532891526049858, 9710696558260003146, ... ]
```


### Nodes

Nodes are stored in a dataframe:

```python
print(graph.nodes)
# node_id                        block_id path_id strand  start    end
# 11484376918084368   6227233701292645975      12   True  87911  88000
# 62802772372552842  14279814672519617104      10   True  30224  30923
# 66596091345916983   1967902255453418588      13   True  27906  28591
# ...                                 ...     ...    ...    ...    ...
# [1042 rows x 5 columns]
```

In this dataframe each node, identified by its `node_id`, is associated with a `block_id` and a `path_id`. The `strand`, `start` and `end` columns indicate the orientation and position of the node in the input path genome.

Each node objects can be retrieved by its id:
```python
node = graph.nodes[11788816242159313242]
print(node)
# block_id    14710008249239879492
# path_id                        0
# strand                      True
# start                       2358
# end                         2552
# Name: 11788816242159313242, dtype: object
```


### Blocks

Blocks encode multiple sequnce alignment. They can be accessed via their id.

```python
block = graph.blocks[14710008249239879492]
print(block)
# block 14710008249239879492, consensus len = 183 bp, n. nodes = 4
```

Each block has a consensus sequence:

```python
print(block.consensus())
# ATATATGGTGCGTTAATTTTTAAACCCTTATTTAATTTC...
```

And an alignment that connects nodes to sequences:

```python
print(block.alignment.generate_alignment())
# { '4174336837421425166': 'ATATATGGTGCGTTAATTTTTAAACCCT...',
#   '8533989107945450583': 'ATATATGGTGCGTTAATTTTTAAACCCT...',
#  '11788816242159313242': 'ATATATGGTGCGTTAATTTTTAAACCCT...',
#  '16194835320646696346': 'ATATATGGTGCGTTAATTTTTAAACCCT...'}
```

More details on alignments are provided in [tutorial 3](t03-block-alignments.md).
