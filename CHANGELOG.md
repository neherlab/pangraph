## 1.1.0

### Make Pangraph CLI logs less verbose

`pangraph build` command now displays only essential information when in "info" verbosity mode (with a single `-v`). To display a more verbose log, add multiple occurences of `-v`: e.g. `-vv`, `-vvv` etc.


### Add progress bar

`pangraph build` command now displays a progress bar which helps to estimate the time required to complete the operation. Add `--no-progress-bar` if you want to hide it.


## 1.0.0

This release introduces several significant updates, including a complete rewrite of the algorithm in Rust to boost both speed and reliability in graph construction. PanGraph can now be compiled as a standalone binary, greatly simplifying installation. The release also contains several **breaking changes**.

### Major changes

- **Graph Construction Algorithm Enhancements**
  - Block merging has been refined and simplified: when merging two blocks, sequences are re-aligned to the block consensus, resulting in more robust alignments.
  - The alignment process is now easily parallelizable, extending beyond the previous limitation of parallel guide tree traversal.
  - Insertions are now placed on the consensus sequence without aligning them to each other, simplifying bookkeeping of alignment variation when merging blocks.

- **Graph JSON Output Format Updates**
  - A new `nodes` dictionary links paths and blocks.
  - Entries in the `paths` dictionary have been simplified to include only lists of nodes.
  - Block `alignments` now use a simpler encoding; insertions are placed on the consensus without alignment, eliminating the need for a separate `gaps` entry.

- **Command Line Interface Modifications**
  - Improved parallelization control with the new `--jobs` flag.
  - The `export` command has been restructured into several subcommands:
    - `export gfa` – Export the graph in GFA format.
    - `export block-consensus` – Export block consensus sequences in a single FASTA file.
    - `export block-sequences` – Export each block’s alignment in separate FASTA files.
    - `export core-genome` – Export the core-genome alignment.
  - The `marginalize` command has been renamed to `simplify`.
  - A new `reconstruct` command has been added to rebuild the input sequences from the graph.
  - The `polish` and `generate` commands have been removed.


## v0.6.0

- added [mmseqs2](https://github.com/soedinglab/MMseqs2) as an alternative alignment kernel that guarantees higher sensitivity at the expense of longer computational time, see [#33](https://github.com/neherlab/pangraph/pull/33).
- updated Docker file to include mmseqs2 in the container.
- updated the documentation, including discussion of alignment kernel sensitivities and examples of application of PanGraph to plasmids by [@liampshaw](https://github.com/neherlab/pangraph/commits?author=liampshaw).
- errors that occur in worker threads are now emitted on the main thread, see [#25](https://github.com/neherlab/pangraph/pull/25).
- fixed bug when using `mash` option see this [commit](https://github.com/neherlab/pangraph/commit/2167c2e9f72b2962ef2e2b9ec1fbe0e16fe0f568)
- fixed a bug in detransitive, see this [commit](https://github.com/neherlab/pangraph/commit/a9651323aba2822d1b1c380a086fae4216c8030d)
- added snakemake pipeline in the `script` folder to perform the analysis published in our [paper](https://github.com/neherlab/pangraph#citing).
- added `-K` option to the `build` command to control kmer length for mmseqs aligner, see this [commit](https://github.com/neherlab/pangraph/commit/0857c36c7c8d11d53e8efab91cf5d18c35685a6e).
- added `fasttree` to docker container and PanX export to docker tests, see [#37](https://github.com/neherlab/pangraph/pull/37).

## v0.5.0

- fix: error with gfa export of fully duplicated paths by @mmolari in [#19](https://github.com/neherlab/pangraph/pull/19)
- GFA export bug fixes by @nnoll in [#28](https://github.com/neherlab/pangraph/pull/28)
- chore: add docker container by @ivan-aksamentov in [#27](https://github.com/neherlab/pangraph/pull/27)
- fix: deal with zero length blocks getting added to segment by @nnoll in [#20](https://github.com/neherlab/pangraph/pull/20)

[Full changelog](https://github.com/neherlab/pangraph/compare/v0.4.1...0.5.0)

## v0.4.1

- Smaller binaries: Artifacts now pulled in as needed.

## v0.4.0

- Marginalize command now (optionally) takes list of strains to project onto.
- Command line arguments and flags can now be mixed in order.
- Export can filter out any duplications.

## v0.3.0

Added command line options:
- Build: `-u` => force sequences to uppercase letters
- Polish: `-c` => preserve case (uses MAFFT command line flag)
Additionally, removed bug associated with sequences mapping as empty intervals to blocks.

Lastly, large improvement to the algorithm's multicore usage by balancing the initial guide tree.

## v0.2.1

Modified CLI to accept input on standard input for all subcommands. This allows for a nicer chaining of pangraph functions from the shell. Additionally, there were many small bugs that are fixed.

## v0.1-alpha

Source code bundled as a relocatable application. Currently only for Linux-based operating systems but intend to release for MacOSX as well.
