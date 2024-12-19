# export command

The `export` command can be used to:
- export in `gfa` v1 format, see [format specifications](https://github.com/GFA-spec/GFA-spec/blob/master/GFA1.md)
- export block consensus sequences to a fasta file
- export block alignments or sequences to multiple fasta files
- export the core-genome alignment

Below is a description for each of the export options.

If `export` is not a great name maybe we could use something like `convert` or `extract`?

## GFA export

The GFA export command can be used to export a pangraph to a GFA file. The command has the following signature:

```bash
# export to gfa
pangraph export gfa graph.json \
    --min-len=100 \
    --min-depth=5 \
    --no-duplicated \
    --include-sequences \
    > graph.gfa
```

It uses the `export_graph_to_gfa` function outlined in `export_gfa.py`.

Optional flags:
- `--min-len` minimum length of a block to be included in the GFA file.
- `--min-depth` minimum depth of a block to be included in the GFA file.
- `--no-duplicated` exclude blocks that are duplicated in any path from the GFA file.
- `--include-sequences` include block sequences in the GFA file.

## export block sequences

The consensus sequence of each block can be exported in a single fasta file with the command:

```bash
pangraph export block-consensus graph.json > block_consensus.fa
```

The command uses the `export_block_consensus` function outlined in `export_block_consensus.py`.

## export block alignments

The export command can also be used to export the sequences for each block.
- `export block-alignments` exports aligned sequences for each block. Note that alignments exclude insertions.
- `export block-sequences` exports the full unaligned sequences for each block, in case the user wants to re-perform alignments with other tools.

These two commands have similar signatures:
```bash
pangraph export block-alignments graph.json -o export_folder
pangraph export block-sequences graph.json -o export_folder
```
The only optional argument is `-o` or `--output` which specifies the output folder. Inside of this folder a series of files `block_<block_id>.fa` are created, each containing the sequences of the block with the corresponding id.

The command used in both is `export_block_sequences` functions outlined in `export_block_sequences.py`, with the flag `aligned` set to `True` or `False` respectively.

## export core-genome alignment

The core-genome alignment can be exported to a fasta file with the command:
```bash
pangraph export core-alignment graph.json \
    --guide-strain=ref_strain \
    --unaligned \
    > core_aln.fa
```

Optional flags are:
- `--guide-strain` specifies the strain to use as a reference for the alignment. Core blocks are ordered and oriented (forward or reverse) according to the reference strain.
- if `--unaligned` is set, then the full core sequences are exported but not aligned, but they should be linearly alignable and can be fed to an external aligner. (This corresponds to setting `aligned=False` in the corresponding python function).

The command uses the `export_core_alignment` function outlined in `export_core_alignment.py`.

## implementation notes

Fasta records are currently implemented [here](https://github.com/neherlab/pangraph/blob/rust/packages/pangraph/src/io/fasta.rs).

Utilities functions contained in `circularize_utils.py` are already implemented [here](https://github.com/neherlab/pangraph/blob/rust/packages/pangraph/src/circularize/circularize_utils.rs)

We can extend the [`apply_edits_to_ref`](https://github.com/neherlab/pangraph/blob/98886771cb20cd4bfe7ce33c52dafc2fc33f6faa/packages/pangraph/src/pangraph/edits.rs#L194) function to output also aligned seqeunces, i.e. sequences without insertions and with gaps at the place of deleted characters.
In the draft code this is done by either adding an `aligned` binary argument to the function. But adding a new function would work as well.