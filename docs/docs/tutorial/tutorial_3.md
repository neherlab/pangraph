---
sidebar_position: 3
---

# Polishing the pangraph and exploring alignments

This third tutorial is focused on block alignments and is divided in three main parts:
- refining the block alignments using the `polish` command.
- exporting the alignment (and related phylogenetic tree) of each block.
- exploring the full set of block alignments using the [panX](https://pangenome.org/)[^1] visualization.

[^1]: Ding, Wei, Franz Baumdicker, and Richard A. Neher. "panX: pan-genome analysis and exploration." Nucleic acids research 46.1 (2018): e5-e5.

!!! note "Requirements"
    Polishing the pangenome graph requires [mafft](https://mafft.cbrc.jp/alignment/software/). Exporting files for the panX visualization requires [fasttree](http://www.microbesonline.org/fasttree). Running the panX visualization requires [node.js](https://nodejs.org/en/).

## Polishing the alignments

The pangenome graph is built by iterative pairwise merges of smaller graphs along a guide tree. At each stage blocks are compared and aligned using only their _consensus_. This shortens considerably the time needed to build a pangraph, but might introduce minor inconsistencies and artifacts in the alignments. These can be removed using the `polish` command which reconstructs the sequences in the block and performs a multiple sequence alignment (e.g. using mafft) (see [Polish](@ref)).

To polish the `ecoli_pangraph.json` file that was produced in the first tutorial we can run the command:

```bash
pangraph polish ecoli_pangraph.json > polished_pangraph.json
```

This should complete in 10-15 minutes on 4 cores. Optionally, a threshold length can be specified with the flag `--length`. In this case only alignments of blocks whose consensus sequence is shorter than this threshold are refined.

Note that the only effect of this command is refining the alignments, but the "topological" structure of the pangraph is left unchanged. As such, polishing a pangraph is only useful when one is directly interested in block alignments.

## Export block alignments and trees

By running the `export` command with the `--export-panX` option (see [Export](@ref)) one can export all the block alignments, together with additional files that are needed for the panX visualization. For our example we can run:

```bash
pangraph export \
    --export-panX \
    --output-directory ecoli_export \
    --no-export-gfa \
    polished_pangraph.json
```

The flag `--no-export-gfa` excludes the export in gfa format that we covered in the first tutorial. This command should run in about 15 minutes. It will generate a folder named `ecoli_export/vis` that contain multiple files. These include:
- the alignment of each block (`geneCluster/*_aln.fa.gz` files)
- the phylogenetic tree generated from each alignment (`geneCluster/*.nwk` files)
- the phylogenetic tree obtained from the concatenated alignment of all core blocks (`strain_tree.nwk`)
- additional files needed for the panX visualization

