# PanGraph Changelog

## v0.7.0

- fasta input files are checked for duplicated records, and white lines between records are tolerated, see [#55](https://github.com/neherlab/pangraph/pull/55).
- PanGraph execution is now deterministic, and same input files always produce the same output, see [#57](https://github.com/neherlab/pangraph/pull/57).
- introduced the `-v` flag in the `build` command. This activates consistency checks at each mergers. In this check we verify that the input genomes can be exactly reconstructed from the graph. See [#57](https://github.com/neherlab/pangraph/pull/57).
- Fixed [#56](https://github.com/neherlab/pangraph/issues/56)

## v0.6.3

- fixed an issue with PanX export by upgrading TreeTools version, see [#52](https://github.com/neherlab/pangraph/issues/52) by @mmolari and @PierreBarrat
- minor upgrades to the analysis scripts following PanGraph's paper reviews, see [#51](https://github.com/neherlab/pangraph/pull/51).

## v0.6.2

- removed the dependency on conda/PyCall and substituted it with [TreeTools](https://github.com/PierreBarrat/TreeTools.jl), by @PierreBarrat and @mmolari, see [#45](https://github.com/neherlab/pangraph/pull/45). This reduces the size of the docker image.
- added `script/config/accnums.json` file with list of accession number for GenBank sequences used for pangraph algorithm validation.
- added `procps` to docker image. This is needed for compatibility of the docker image with nextflow. By @plaquette, see [#48](https://github.com/neherlab/pangraph/pull/48).
- minor improvements in the analysis pipeline, see [#49](https://github.com/neherlab/pangraph/pull/49).

## v0.6.1

- added `pangraph version` command that prints PanGraph's version on stderr, by @mmolari and @ivan-aksamentov, see [#40](https://github.com/neherlab/pangraph/pull/40).
- fix: wrong PanGraph's package version tag in `Project.toml`.
- added `tools/release.sh` script to automate the release process, by @mmolari and @ivan-aksamentov, see [#41](https://github.com/neherlab/pangraph/pull/41).

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
