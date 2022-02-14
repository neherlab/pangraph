# Building a pangraph

This short tutorial will walk you through the process of generating a pangraph from a set of bacterial genomes. We will also cover how to export the generated pangraph file into other formats.

Simply put, a **pangenome graph** (or _pangraph_ for short) is a compressed representation of a set of genomes, in which alignable regions are saved as _blocks_ (or _pancontigs_) and genomes are represented as _paths_, i.e. list of blocks. 

![img](./../assets/pangraph_scheme.png)


## Requirements

The tutorial requires you to have the `pangraph` command available in your path. Instructions on how to install pangraph can be found in [Installation](@ref).

For this tutorial we will use a small dataset containing full chromosomes of 10 _Escherichia Coli_ strains (source: GenBank). For convenience this dataset is available in the pangraph repository (`example_dataset/ecoli.fa.gz`), and can be downloaded with the command:

```bash
wget https://github.com/nnoll/pangraph/raw/feat/example-dataset/example_datasets/ecoli.fa.gz
```

This is a single fasta file containing 10 entries. We chose to include no plasmids in the example dataset. Notice that it is not necessary for all of the data to be packed in a single fasta file. One can also pass multiple fasta files (optionally gzipped) as input to the command to build a pangraph.

## Building the pangraph

As a first step, we will build a pangraph object from the DNA of the 10 chromosomes. This can be done using the command `build` (see [Build](@ref)):

```bash
pangraph build --circular ecoli.fa.gz > ecoli_pangraph.json
```

On a consumer laptop the command should complete in around 10 minutes. The option `--circular` is used when passing circular DNA sequences, like the bacterial chromosomes that we consider here.

The result is a `ecoli_pangraph.json` file that contains two main entries: `paths` and `blocks`. As represented in the image above, blocks contain information on the nucleotide sequence, while paths are compressed representation for genomes as lists of blocks. Importantly, each block is assigned an unique random id composed of 10 capital letters. Below is a simplified view of the structure of the `ecoli_pangraph.json` file.

```json
{
    "paths": [
        {
            "name": "NZ_CP010242",
            "blocks": [ { "id": "NFTNKNMFIC", ... }, { "id": "YTLSRRNNGL", ... }, ... ],
            ...
        },
        {
            "name": "NC_009800",
            "blocks": [  { "id": "AYYUXVXZXB", ... },  { "id": "YTLSRRNNGL", ... }, ... ],
            ...
        }
        ...
    ],
    "blocks": [
        { "id": "KZJIDOXBAV", "sequence": "AAGGTGGGTAATCATTTTGATAAGTGAT...", ... },
        { "id": "UOFDTEUSWC", "sequence": "GTTTTAATGCCAGCAAAAATGGTGAATT...", ... },
        ...
    ]
}
```

Each entry in `path` has two main properties: the `name`, corresponding to the sequence identifier in the input fasta file, and the `blocks` list. The latter is a representation of the genome as a list of blocks, each one identified by its unique `id`.

Each entry in the `blocks` lists corresponds instead to a different block. The two main properties of each entry are the block `id` and the consensus `sequence` of the block.

More details on the structure of this `json` file will be covered in the next tutorial section.


## Exporting the pangraph

The pangraph object can also be exported in other more common formats using the command `export` (see [Export](@ref)).

```bash
pangraph export \
    --no-duplications \
    --output-directory ecoli_export \
    ecoli_pangraph.json
```

This will create a folder named `ecoli_export` that contains two files.

- `pangraph.fa`: a fasta file containing the consensus sequence for each block.
- `pangraph.gfa`: a [Graphical Fragment Assembly](https://github.com/GFA-spec/GFA-spec) file that contains a representation of the pangenome graph structure.

The latter can be visualized using [Bandage](https://rrwick.github.io/Bandage/). The option `--no-duplications` causes the export function to avoid including duplicated blocks in the graph representation (they are instead exported as isolated blocks). In our experience this results in a less "tangled" visual representation. Below is how the Bandage visualization of this example pangraph looks like. Blocks are colored by frequency, with common blocks (appearing in many different chromosomes) in red and rare blocks (appearing in only a few chromosomes) in black. 

![img](./../assets/bandage_ecoli_full.png)
