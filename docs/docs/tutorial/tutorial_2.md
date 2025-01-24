---
sidebar_position: 2
---

# The structure of Pangraph output file

In this second part of the tutorial we will explore in more detail the content of the `json` output file produced by the [`build` command](../reference#pangraph-build).
As an example, we will use snippets from the `graph.json` file that was produced in the previous section of the tutorial.


## The structure of `graph.json`

As discussed in the [previous tutorial section](./tutorial_1.md#what-is-a-pangraph), the three main entries of pangraph output file are `paths`, `blocks` and `nodes`.

- each entry in the `paths` list encodes one of the nucleotide sequences that were given as input to the `build` command, represented as a list of nodes (i.e. particular instances of a block)
- each entry in the `blocks` list represents an alignable set of homologous sequences. A block contains the consensus of all of these sequences, together with information to reconstruct the full alignment. Each entry in the alignment is represented by a `node`.
- the `nodes` list represents the connection between blocks and paths. Each node entry contains information on which block and path the node is assigned to.

We will explore each of these categories separately.


## Paths

A path object has the following structure:

```json
{
    "id": 1,
    "name": "NZ_CP010150",
    "nodes": [10429785587629589393, 10765941013351965021, 7771937209474314297, ...],
    "circular": true,
    "tot_len": 4827779
},
```

The two main properties of a path are `name` and `nodes`. The `name` of the path corresponds to the sequence identifier in the input fasta file.
`nodes` contains the ordered list of node ids that make up the path.

Here is a complete list containing a description of every entry in the path object:

- `id`: numerical id of the path. Each path is assigned a unique progressive id when building the graph.
- `name` : the path sequence identifier, as specified in the input fasta file.
- `nodes` : the ordered list of node ids that make up the path.
- `circular` : indicates whether the considered sequence is circular (e.g. plasmid) or not. This is controlled by the `--circular` flag of the `build` command.
- `tot_len` : the total length of the path, in nucleotides.


## Nodes

Nodes constitue the connection between block and paths. A node indicartes a particular occurrence of a block in a path sequence. A node object has the following structure:

```json
{
    "id": 10429785587629589393,
    "block_id": 9245376340613946,
    "path_id": 1,
    "strand": "+",
    "position": [356656, 359732]
},
```

The properties of a node are:
- `id` : a unique random numerical id assigned to the node.
- `block_id` : the id of the block that the node belongs to.
- `path_id` : the id of the path that the node is part of.
- `strand` : indicates whether on the original input sequence the consensus sequence of the block appears in the forward (`+`) or reverse (`-`) strand.
- `position` : start and end position of the node on the input sequence. Positions are always in 0-based numbering with right extreme excluded. They are also based on the forward strand (with beginning < end). The only exception is when a block wraps around the end of a circular sequence. In this case the node start position (close to the end of the genome) is higher than the node end position (close to the beginning of the genome).


## Blocks

Blocks encode alignments of homologous sequences across the input genomes.
Each block object has the following structure:

```json
{
    "id": 9783460543474296593,
    "consensus": "AACGGCAATATCTGCCACAAA...",
    "alignments": {
        "5252840658835653895": {
            "subs": [ ... ],
            "dels": [ ... ],
            "inss": [ ... ]
        },
        "7608024242617339186": {
            "subs": [ ... ],
            "dels": [ ... ],
            "inss": [ ... ]
        },
        ...
    }
}
```

Each block contains the following properties:
- `id` : a unique random numerical id assigned to the block.
- `consensus` : the consensus sequence of the block.
- `alignments` : a dictionary that contains information to reconstruct the sequence alignment. Keys are node ids, while values are a list of variations (substitutions, deletions, insertions) that need to be applied to the consensus to obtain the node sequence.


### How alignments are encoded

Each object in the alignment dictionary contains three entries: `subs`, `dels`, and `inss`, encoding respectively substitutions, deletions, and insertions. Here is an example of an entry in the `alignments` dictionary:

```json
"5252840658835653895" : {
    "subs": [
        {
            "pos": 15,
            "alt": "A"
        },
        {
            "pos": 35,
            "alt": "G"
        }
    ],
    "dels": [
        {
            "pos": 10,
            "len": 4
        }
    ],
    "inss": [
        {
            "pos": 3,
            "seq": "CTT"
        }
    ]
}
```

- `subs` contains a list of substitutions. Each substituion contains a `pos` field indicating the position of the substitution in the consensus sequence, and an `alt` field indicating the alternative nucleotide.
- `dels` contains a list of deletions. Each deletion contains a `pos` field indicating the starting position of the deletion in the consensus sequence (0-based), and a `len` field indicating the length of the deletion.
- `inss` contains a list of insertions. Each insertion contains a `pos` field, indicating the nucleotide _before_ which the insertion should be placed, and a `seq` field, containing the sequence to be inserted.

Below is a schematic representation of how these variations are applied to the consensus sequence of a block to obtain the sequence of a node.

![img](./../assets/t2_alignment_reconstruction.png)

As discussed in the [next section](./tutorial_3.md), the different sequences of a block can be reconstructed in two ways:
- as **node sequences**. In this case sequences are not aligned, but each entry corresponds to the exact sequence of a node, with all variations applied.
- as a **multiple sequence alignment**. In this case sequences are aligned, but insertions are omitted. 


## A look at the length and frequency of blocks

Having completed this part of the tutorial, it is now possible to access directly the rich information contained in the pangraph output format.

As a simple example, we take the `ecoli_pangraph.json` file and extract from each block the length of its consensus sequence, and the number of unique E.Coli strains that contain the block. We then plot the cumulative distribution of block lengths, and the frequency distribution of blocks (weighted by length).

![img](../assets/ecoli_blockstats.png)

The blocks present in this pangraph have widely varying size, with some blocks being only some hundreds of nucleotides long, and others spanning tens of kbps.[^1]

[^1]: The almost complete absence of blocks shorter than 100 bps is due to the default value of the `--len` option in the `build` command, which is set to 100. See [`build` command](../reference#pangraph-build) for more information.

The block frequency distribution shows a typical bimodal pattern, with an abundance of "core" blocks (blocks that are present in all the 10 considered chromosomes, cumulative length of more than 3.5 Mbps) and rare blocks present in only one strain (cumulative length of more than 2.5 Mbps).
This graph can be reproduced by this python snippet:
```python
import json
import matplotlib.pyplot as plt

with open('ecoli_pangraph.json') as fh:
    G = json.load(fh)

blocks = [{'length':len(x['sequence']), 'depth':len(x['positions'])} for x in G['blocks']]
plt.hist([x['depth'] for x in blocks], weights=[x['length'] for x in blocks], bins=range(1,11))
```

