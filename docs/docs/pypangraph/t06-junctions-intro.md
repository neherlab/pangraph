---
sidebar_position: 8
---

# Junctions: accessory genome structure

When comparing closely related bacterial genomes, the core genome is largely syntenic — long stretches of conserved blocks appear in the same order across isolates. Between these conserved blocks, segments of **accessory** DNA (insertions, deletions, mobile elements) vary from one genome to another. A **junction** captures one such segment together with its flanking core blocks, providing a natural unit for comparing structural variation across isolates.

In this tutorial we introduce the junction concept and the `BackboneJunctions` class, which is the main entry point for junction analysis in pypangraph.

## What is a junction?

To define junctions, we first identify **backbone blocks**: core blocks (present exactly once in every genome) whose consensus length is at least `L_thr` base pairs (default 500). These long, universally conserved blocks form a stable "spine" through the pangenome.

A **junction** is the segment on a given genome's path between two consecutive backbone blocks:

```
   [left backbone block] ── accessory center ── [right backbone block]
         (core)           (0 or more blocks)          (core)
```

The two flanking backbone blocks define an **edge**, identified by a string encoding the block IDs and strand of the flanks (e.g. `"124231456905500231_f__1603316146112203317_f"`). Different isolates that share the same edge can have different accessory content in the center — different blocks, different lengths, or no accessory blocks at all.

:::info edge string format

Edge strings have the format `<left_block_id>_<strand>__<right_block_id>_<strand>`, where `_f` means forward strand and `_r` means reverse. Each edge has a single canonical ID: the lexicographically smaller of the two possible orientations (the edge vs. its reverse complement).

:::

## Creating BackboneJunctions

We start by loading a graph and creating a `BackboneJunctions` object:

```python
import pypangraph as pp

graph = pp.Pangraph.from_json("plasmids.json")
print(graph)
# pangraph object with 15 paths, 137 blocks and 1042 nodes

bj = pp.junctions.BackboneJunctions(graph, L_thr=500)
```

The `L_thr` parameter controls which core blocks count as backbone. Only core blocks with consensus length >= `L_thr` are used as junction boundaries. Shorter core blocks are absorbed into the accessory center of junctions. This avoids splitting junctions at very small core blocks that may not represent stable genomic landmarks.

Results are lazily computed: the first call to any method triggers the path splitting, after which results are cached.

## Listing edges and junctions

```python
edges = bj.edges()
print(f"Number of edges: {len(edges)}")
# Number of edges: 20

print(edges[:3])
# ['124231456905500231_r__865151745502309237_r',
#  '124231456905500231_f__1603316146112203317_f',
#  '1603316146112203317_f__8434022508348362741_f']
```

We can also list junctions for a specific isolate:

```python
juncs = bj.junctions_for("RCS48_p1")
print(f"Number of junctions: {len(juncs)}")
# Number of junctions: 20

# inspect one junction
j = juncs[0]
print(j)
# [865151745502309237|+|n...] <-- [14710...|+|n...]_[8000...|+|n...]_... --> [124231...|+|n...]
```

Each junction has a `left` flank, a `center` (a path of accessory blocks), and a `right` flank.

## Junction dataframe

The `dataframe()` method returns a pivot table summarizing accessory content per isolate and edge, together with per-edge statistics:

```python
jdf, stats_df = bj.dataframe()
print(jdf.shape)
# (15, 20)

print(jdf.iloc[:3, :3])
# edge       124231..._r__865151..._r  124231..._f__160331..._f  160331..._f__843402..._f
# iso
# RCS100_p1                    29083.0                       0.0                     327.0
# RCS29_p1                     14961.0                       0.0                     329.0
# RCS33_p1                     14961.0                       0.0                     329.0
```

Rows are isolates, columns are edges, and values are the total accessory block length in each junction. `NaN` means the isolate does not have that edge (which would indicate a rearrangement). The `stats_df` dataframe contains per-edge statistics such as frequency and diversity — see [junction statistics](t07-junction-stats.md) for details.

## Co-orientation

Junctions for the same edge may appear in opposite orientations on different genomes: one genome reads the junction left-to-right, while another traverses it right-to-left. The `BackboneJunctions` class handles co-orientation automatically when comparing junctions or extracting sequences. This is transparent to the user in most cases.

<details>
<summary>How co-orientation works</summary>

The canonical edge string defines a reference direction. When a junction's left flank matches the canonical edge's left block (same block ID and strand), the junction is in canonical orientation. Otherwise, it is inverted.

For operations like computing path categories (to determine if two isolates have the same accessory content) or extracting sequences, junctions are co-oriented to the canonical direction before comparison. This ensures that the same accessory content is recognized as identical regardless of which direction a genome traverses it.

</details>

## Next steps

In the following tutorials we explore junction analysis in more detail:

- [Junction statistics](t07-junction-stats.md) — characterize structural variation across isolates
- [Junction genomic positions](t08-junction-positions.md) — map junctions to genomic coordinates
- [Extracting junction sequences](t09-junction-sequences.md) — export co-oriented sequences for downstream analysis
