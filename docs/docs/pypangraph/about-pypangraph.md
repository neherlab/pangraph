---
sidebar_position: 1
---

# About PyPangraph

PyPanGraph is a Python package to facilitate exploration and analysis of [PanGraph](https://github.com/neherlab/pangraph) output JSON files.

PyPangraph can be installed following [these instructions](installation.md).

Below you'll find some simple usage of PyPangraph. For a more complete guide you can follow the [tutorials](t01-load-graph.md).

```python
import pypangraph as pp

# load a graph
graph = pp.Pangraph.from_json("graph.json")
# pangraph object with 15 paths, 137 blocks and 1042 nodes

# recover a specific path with its identifier
path = graph.paths["RCS48_p1"]
# path object | name = RCS48_p1, n. nodes = 60, length = 80596 bp

# extract a block alignment
block = graph.blocks[124231456905500231]
# block 124231456905500231, consensus len = 183 bp, n. nodes = 4

aln = block.to_biopython_alignment()
# Alignment with 15 rows and 2932 columns
# TTCTGCAATTGAGTCTTGTATGCCCCCATAACAGCACTAAATAA...
# TTCTGTAATTGAGTCTTGTATGCCCCCATAACAGCACTAAATAA...
# TTCTGTAATTGAGTCTTGTATGCCCCCATAACAGCACTAAATAA...
# TTCTGCAATTGAGTCTTGTATGCCCCCATAACAGCACTAAATAA...
# ...

# get blocks statistics (length, copy number...)
stats_df = graph.to_blockstats_df()
# block_id              count  n_strains  duplicated   core   len
# 124231456905500231       15         15       False   True  2202
# 149501466629434994        2          2       False  False   210
# 279570196774736738        2          2       False  False  1308
# ...                     ...        ...         ...    ...   ...
```