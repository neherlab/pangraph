---
sidebar_position: 9999
---
# Reference


This document contains the automatically generated reference documentation for command-line arguments of the latest version of Pangraph CLI.

If you have Pangraph CLI installed, you can type `pangraph --help` to read the latest documentation for your installed version of Pangraph. To generate this document in markdown format, run `pangraph help-markdown > reference.md`
  

**Command Overview:**

* [`pangraph`↴](#pangraph)
* [`pangraph build`↴](#pangraph-build)
* [`pangraph export`↴](#pangraph-export)
* [`pangraph export gfa`↴](#pangraph-export-gfa)
* [`pangraph export block-consensus`↴](#pangraph-export-block-consensus)
* [`pangraph export block-sequences`↴](#pangraph-export-block-sequences)
* [`pangraph export core-genome`↴](#pangraph-export-core-genome)
* [`pangraph simplify`↴](#pangraph-simplify)
* [`pangraph reconstruct`↴](#pangraph-reconstruct)
* [`pangraph schema`↴](#pangraph-schema)
* [`pangraph completions`↴](#pangraph-completions)
* [`pangraph help-markdown`↴](#pangraph-help-markdown)

## `pangraph`

Bioinformatic toolkit to align large sets of closely related genomes into a graph data structure.

Finds homology amongst large collections of closely related genomes. The core of the algorithm partitions each genome into pancontigs that represent a sequence interval related by vertical descent. Each genome is then an ordered walk along pancontigs; the collection of all genomes form a graph that captures all observed structural diversity. The tool useful to parsimoniously infer horizontal gene transfer events within a community; perform comparative studies of genome gain, loss, and rearrangement dynamics; or simply to compress many related genomes.


Publication: "PanGraph: scalable bacterial pan-genome graph construction. Nicholas Noll, Marco Molari, Richard Neher. bioRxiv 2022.02.24.481757; doi: https://doi.org/10.1101/2022.02.24.481757"

Documentation: https://pangraph.readthedocs.io/en/stable/

Source code:https://github.com/neherlab/pangraph

Questions, ideas, bug reports: https://github.com/neherlab/pangraph/issues

**Usage:** `pangraph [OPTIONS] <COMMAND>`

###### **Subcommands:**

* `build` — Align genomes into a multiple sequence alignment graph
* `export` — Export a pangraph to a chosen file format(s)
* `simplify` — Compute all pairwise marginalizations of a multiple sequence alignment graph
* `reconstruct` — Reconstruct input fasta sequences from graph
* `schema` — Generate JSON schema for Pangraph file format
* `completions` — Generate shell completions
* `help-markdown` — Print command-line reference documentation in Markdown format

###### **Options:**

* `--verbosity <VERBOSITY>` — Set verbosity level of console output

  Default value: `warn`

  Possible values: `off`, `error`, `warn`, `info`, `debug`, `trace`

* `--silent` — Disable all console output. Same as `--verbosity=off`
* `-v`, `--verbose` — Make console output more verbose. Add multiple occurrences to increase verbosity further
* `-q`, `--quiet` — Make console output more quiet. Add multiple occurrences to make output even more quiet
* `-j`, `--jobs <JOBS>` — Number of processing jobs. If not specified, all available CPU threads will be used



## `pangraph build`

Align genomes into a multiple sequence alignment graph

**Usage:** `pangraph build [OPTIONS] [INPUT_FASTAS]...`

###### **Arguments:**

* `<INPUT_FASTAS>` — Path(s) to zero, one or multiple FASTA files with input sequences. Multiple records within one file are treated as separate genomes.

   Accepts plain or compressed FASTA files. If a compressed fasta file is provided, it will be transparently decompressed. Supported compression formats: `gz`, `bz2`, `xz`, `zstd`. Decompressor is chosen based on file extension. If there's multiple input files, then different files can have different compression formats.

   If no input files provided, the plain fasta input is read from standard input (stdin).

   See: https://en.wikipedia.org/wiki/FASTA_format

###### **Options:**

* `-o`, `--output-json <OUTPUT_JSON>` — Path to output JSON file with resulting pangraph.

   If the provided file path ends with one of the supported extensions: "gz", "bz2", "xz", "zst", then the file will be written compressed. If the required directory tree does not exist, it will be created.

   Use "-" to write the uncompressed data to standard output (stdout). This is the default, if the argument is not provided.

  Default value: `-`
* `-l`, `--len <INDEL_LEN_THRESHOLD>` — Minimum block size for alignment graph (in nucleotides)

  Default value: `100`
* `-a`, `--alpha <ALPHA>` — Energy cost for splitting a block during alignment merger. Controls graph fragmentation, see documentation

  Default value: `100`
* `-b`, `--beta <BETA>` — Energy cost for diversity in the alignment. A high value prevents merging of distantly-related sequences in the same block, see documentation

  Default value: `10`
* `-s`, `--sensitivity <SENSITIVITY>` — Used to set pairwise alignment sensitivity for minimap aligner. Corresponds to option -x asm5/asm10/asm20 in minimap2

  Default value: `10`
* `-K`, `--kmer-length <KMER_LENGTH>` — Sets kmer length for mmseqs2 aligner
* `-c`, `--circular` — Toggle if input genomes are circular
* `-x`, `--max-self-map <MAX_SELF_MAP>` — Maximum number of alignment rounds to consider per pairwise graph merger

  Default value: `100`
* `-k`, `--alignment-kernel <ALIGNMENT_KERNEL>` — Backend to use for pairwise genome alignment

  Default value: `minimap2-lib`

  Possible values: `minimap2-lib`, `minimap2-cli`, `mmseqs`

* `-f`, `--verify` — Sanity check: after construction verifies that the original sequences can be reconstructed exactly from the resulting pangraph. Raises an error otherwise



## `pangraph export`

Export a pangraph to a chosen file format(s)

**Usage:** `pangraph export <COMMAND>`

###### **Subcommands:**

* `gfa` — Export to GFA v1 format
* `block-consensus` — Export block consensus sequences to a fasta file
* `block-sequences` — Export aligned or unaligned sequences for each block. Note that alignments exclude insertions
* `core-genome` — Export the core-genome alignment



## `pangraph export gfa`

Export to GFA v1 format

See GFA v1 format specifications: https://github.com/GFA-spec/GFA-spec/blob/master/GFA1.md

**Usage:** `pangraph export gfa [OPTIONS] [INPUT_JSON]`

###### **Arguments:**

* `<INPUT_JSON>` — Path to a pangraph file (native json).

   Accepts plain or compressed files. If a compressed file is provided, it will be transparently decompressed. Supported compression formats: `gz`, `bz2`, `xz`, `zstd`. Decompressor is chosen based on file extension.

   If no path provided, the uncompressed input is read from standard input (stdin).

###### **Options:**

* `-o`, `--output <OUTPUT>` — Path to output GFA file.

   See GFA v1 format specifications: https://github.com/GFA-spec/GFA-spec/blob/master/GFA1.md

   If the provided file path ends with one of the supported extensions: "gz", "bz2", "xz", "zst", then the file will be written compressed. If the required directory tree does not exist, it will be created.

   Use "-" to write the uncompressed to standard output (stdout). This is the default, if the argument is not provided.

  Default value: `-`
* `--minimum-length <MINIMUM_LENGTH>` — Blocks below this length cutoff will not be exported
* `--maximum-length <MAXIMUM_LENGTH>` — Blocks above this length cutoff will not be exported
* `--minimum-depth <MINIMUM_DEPTH>` — Blocks below this depth cutoff will not be exported
* `--maximum-depth <MAXIMUM_DEPTH>` — Blocks above this depth cutoff will not be exported
* `--include-sequences` — Include block sequences in the GFA file
* `--no-duplicated` — Exclude blocks that are duplicated in any path



## `pangraph export block-consensus`

Export block consensus sequences to a fasta file

**Usage:** `pangraph export block-consensus [OPTIONS] [INPUT_JSON]`

###### **Arguments:**

* `<INPUT_JSON>` — Path to a pangraph file (native json).

   Accepts plain or compressed files. If a compressed file is provided, it will be transparently decompressed. Supported compression formats: `gz`, `bz2`, `xz`, `zstd`. Decompressor is chosen based on file extension.

   If no path provided, the uncompressed input is read from standard input (stdin).

###### **Options:**

* `-o`, `--output <OUTPUT>` — Path to output FASTA file.

   See: https://en.wikipedia.org/wiki/FASTA_format

   If the provided file path ends with one of the supported extensions: "gz", "bz2", "xz", "zst", then the file will be written compressed. If the required directory tree does not exist, it will be created.

   Use "-" to write the uncompressed to standard output (stdout). This is the default, if the argument is not provided.

  Default value: `-`



## `pangraph export block-sequences`

Export aligned or unaligned sequences for each block. Note that alignments exclude insertions

**Usage:** `pangraph export block-sequences [OPTIONS] --output <OUTPUT> [INPUT_JSON]`

###### **Arguments:**

* `<INPUT_JSON>` — Path to a pangraph file (native json).

   Accepts plain or compressed files. If a compressed file is provided, it will be transparently decompressed. Supported compression formats: `gz`, `bz2`, `xz`, `zstd`. Decompressor is chosen based on file extension.

   If no path provided, the uncompressed input is read from standard input (stdin).

###### **Options:**

* `-o`, `--output <OUTPUT>` — Path to directory to write output FASTA files to

   See: https://en.wikipedia.org/wiki/FASTA_format
* `--unaligned` — If set, then the full block sequences are exported but not aligned



## `pangraph export core-genome`

Export the core-genome alignment

**Usage:** `pangraph export core-genome [OPTIONS] --guide-strain <GUIDE_STRAIN> [INPUT_JSON]`

###### **Arguments:**

* `<INPUT_JSON>` — Path to a pangraph file (native json).

   Accepts plain or compressed files. If a compressed file is provided, it will be transparently decompressed. Supported compression formats: `gz`, `bz2`, `xz`, `zstd`. Decompressor is chosen based on file extension.

   If no path provided, the uncompressed input is read from standard input (stdin).

###### **Options:**

* `-o`, `--output <OUTPUT>` — Path to output FASTA file.

   See: https://en.wikipedia.org/wiki/FASTA_format

   If the provided file path ends with one of the supported extensions: "gz", "bz2", "xz", "zst", then the file will be written compressed. If the required directory tree does not exist, it will be created.

   Use "-" to write the uncompressed to standard output (stdout). This is the default, if the argument is not provided.

  Default value: `-`
* `--guide-strain <GUIDE_STRAIN>` — Specify the strain to use as a reference for the alignment. Core blocks are ordered and oriented (forward or reverse) according to the reference strain
* `--unaligned` — If set, then the full core sequences are exported but not aligned.

   They should be linearly alignable and can be fed to an external aligner.



## `pangraph simplify`

Compute all pairwise marginalizations of a multiple sequence alignment graph

**Usage:** `pangraph simplify [OPTIONS] [INPUT]`

###### **Arguments:**

* `<INPUT>` — Path to Pangraph JSON.

   Accepts plain or compressed file. If a compressed file is provided, it will be transparently decompressed. Supported compression formats: `gz`, `bz2`, `xz`, `zstd`. Decompressor is chosen based on file extension.

   If no input file provided, the uncompressed input is read from standard input (stdin).

###### **Options:**

* `-o`, `--output <OUTPUT>`

  Default value: `-`
* `-s`, `--strains <STRAINS>` — Isolates to project onto: collapse the graph to only blocks contained by paths of the given isolates. List of strain names, comma-delimited without spaces



## `pangraph reconstruct`

Reconstruct input fasta sequences from graph

**Usage:** `pangraph reconstruct [OPTIONS] [INPUT_GRAPH]`

###### **Arguments:**

* `<INPUT_GRAPH>` — Path to a pangenome graph file in JSON format.

   Accepts plain or compressed FASTA files. If a compressed fasta file is provided, it will be transparently decompressed. Supported compression formats: `gz`, `bz2`, `xz`, `zstd`. Decompressor is chosen based on file extension. If there's multiple input files, then different files can have different compression formats.

   If no input files provided, the plain fasta input is read from standard input (stdin).

###### **Options:**

* `-o`, `--output-fasta <OUTPUT_FASTA>` — Path to output FASTA file with reconstructed sequences.

   If the provided file path ends with one of the supported extensions: "gz", "bz2", "xz", "zst", then the file will be written compressed. If the required directory tree does not exist, it will be created.

   Use "-" to write the uncompressed to standard output (stdout). This is the default, if the argument is not provided. See: https://en.wikipedia.org/wiki/FASTA_format

  Default value: `-`
* `-f`, `--verify <VERIFY>` — Path to the FASTA file with sequences to check the reconstructed sequences against. If this argument is provided, then the sequences are not being printed to standard output (stdout) as usual. Instead, if any differences are detected, a diff will be printed between the expected (original) sequence and reconstructed sequence.

   Accepts plain or compressed FASTA files. If a compressed fasta file is provided, it will be transparently decompressed. Supported compression formats: `gz`, `bz2`, `xz`, `zstd`. Decompressor is chosen based on file extension. If there's multiple input files, then different files can have different compression formats.

   Use "-" to read uncompressed FASTA from standard input (stdin).



## `pangraph schema`

Generate JSON schema for Pangraph file format

**Usage:** `pangraph schema [OPTIONS]`

###### **Options:**

* `-o`, `--output-pangraph-schema <OUTPUT_PANGRAPH_SCHEMA>` — Path to output JSON or YAML file to write generated JSON Schema definitions to.

   Use "-" to write the uncompressed data to standard output. This is the default, if this argument is not provided.

   See: https://json-schema.org

  Default value: `-`



## `pangraph completions`

Generate shell completions.

This will print the completions file contents to the console. Refer to your shell's documentation on how to install the completions.

Example for Ubuntu Linux:

pangraph completions bash > ~/.local/share/bash-completion/pangraph

**Usage:** `pangraph completions [SHELL]`

###### **Arguments:**

* `<SHELL>` — Name of the shell to generate appropriate completions

  Default value: `bash`

  Possible values: `bash`, `elvish`, `fish`, `fig`, `powershell`, `zsh`




## `pangraph help-markdown`

Print command-line reference documentation in Markdown format

**Usage:** `pangraph help-markdown`



<hr/>

<small><i>
    This document was generated automatically by
    <a href="https://crates.io/crates/clap-markdown"><code>clap-markdown</code></a>.
</i></small>

