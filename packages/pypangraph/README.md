# PyPangraph

This repository contains a collection of utilities to load, explore and analyze pangrenome graphs produced by [PanGraph](https://github.com/neherlab/pangraph).

The package can be installed via pip, see [the documentation](https://docs.pangraph.org/pypangraph/installation):
```bash
pip install pypangraph
```

Below are some examples showcasing some of the main functions in the package. More detailed information and examples can be found in the [documentation](https://docs.pangraph.org/category/pypangraph).

## Loading and interacting with pangraph objects:

A pangraph object can be loaded from a json file with:

```python
# load the library
import pypangraph as pp
# this returns a pangraph object
graph = pp.Pangraph.from_json("path/to/pangraph.json")
```

Pangraph objects have three main properties: `blocks`, `paths`, and `nodes`.

- `blocks` encode for alignments of homologous sequences across genomes. Each entry in the alignment is a `node`.
- `paths` encode genomes as a list of `nodes`.
- `nodes` connect paths and blocks.

See the documentation for more details on this data structure.

These elements are contained in different properties of the graph:

```python
graph.blocks # dictionary of block ids -> block objects
graph.paths # dictionary of path ids -> path objects
graph.nodes # dataframe of nodes
```

And can be accessed by either iterating through the items or via their ids.


### Block object

Block object have these properties/methods:
```python
# get a block by its id
block = graph.blocks[12252014572476775186]

block.id # unique id of the block
block.consensus() # block consensus sequence
block.depth() # n. of nodes in the block alignment
block.to_biopython_alignment() # returns a biopython alignment object
```

### Path object


```python
# paths can be accessed by their name
path = graph.paths['NZ_CP014647']
# or by their numerical id
path = graph.paths.list[4]

path.id # path numerical id
path.name # path name
path.nuc_len # total length of the path in nucleotides
path.circular # whether the path is circular or linear
path.nodes # list of node ids in the path
```

### Node object

The `nodes` property of the pangraph is a pandas dataframe with the following columns:

```python
graph.nodes
                    
# node_id                           block_id path_id strand  start    end
# 11484376918084368      6227233701292645975      12   True  87911  88000
# 31660532043830364     12252014572476775186       7  False  88597  91675
# 35440216894469496      5326177636996110751      12   True  12629  13292
# ...
```

- `node_id` is the unique id of the node
- `block_id` is the id of the block the node belongs to
- `path_id` is the id of the path the node belongs to
- `strand` is a boolean indicating whether the node occurs in forward (`True`) or reverse (`False`) orientation in the path.
- `start` and `end` are the start and end positions of the node in the input genome.

Nodes can be accessed by their id:

```python
node = graph.nodes[12252014572476775186]
```

### Block statistics

The function `to_blockstats_df()` returns a pandas dataframe containing summary statistics on the blocks, indexed by `block id`:

```python
graph.to_blockstats_df()

# block_id              count  n_strains  duplicated   core   len
# 124231456905500231       15         15       False   True  2202
# 149501466629434994        2          2       False  False   210
# 279570196774736738        4          2        True  False  1308
# ...                     ...        ...         ...    ...   ...
```

- `count` is the number of entries in the block alignment, i.e. the total number of times the block appears in all paths.
- `n_strains` is the number of unique paths in which the block appears. It can be different from `count` if the block is duplicated in some paths.
- `duplicated` indicates whether the block is duplicated in any paths.
- `core` indicates whether the block is core, i.e. present exactly once in every path.
- `len` is the length of the block consensus sequence in basepairs.

### Block count matrix


The function `to_blockcount_df()` returns a pandas dataframe whose columns are path ids, and indices are block ids. The entries are the number of times a block is present in a given strain.

```python
graph.to_blockcount_df()
# path_id               RCS48_p1  RCS49_p1  RCS64_p2  ...
# block_id
# 124231456905500231           1         1         1  ...
# 149501466629434994           0         1         0  ...
# 279570196774736738           0         0         2  ...
# ...                        ...       ...       ...  ...
```