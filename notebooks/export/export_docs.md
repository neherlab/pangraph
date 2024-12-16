# export command

I imagine that the `export` command might be used to export in different formats. For example:
```bash
# export to gfa
pangraph export gfa graph.json -o export_folder

# export block consensus sequences
pangraph export block-consensus graph.json > block_consensus.fa

# export unaligned sequences, to be used in other tools
pangraph export unaligned-sequences graph.json -o export_folder

# export alignments for each block
pangraph export block-alignments graph.json -o export_folder

# export core-genome alignment
pangraph export core-alignment graph.json \
    --guide-strain=ref_strain \
    > core_aln.fa
```

Some design decisions that we need to take:
- we do not align insertions. Exported alignments do not have insertions.
- what are the id of fasta records for block alignments? Node id (and strain name)?
- what is the order of blocks in the core-genome alignment? Let the use specify a `guide-strain`.

## blocks to fasta records

Fasta records are currently implemented [here](https://github.com/neherlab/pangraph/blob/rust/packages/pangraph/src/io/fasta.rs).

We can extend the [`apply_edits_to_ref`](https://github.com/neherlab/pangraph/blob/98886771cb20cd4bfe7ce33c52dafc2fc33f6faa/packages/pangraph/src/pangraph/edits.rs#L194) function to output also aligned seqeunces, i.e. sequences without insertions and with gaps at the place of deleted characters.
In the draft code this is done by either adding an `aligned` binary argument to the function. But adding a new function would work as well.