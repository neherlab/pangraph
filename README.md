# pangraph

![readthedocs](https://nnoll.github.io/pangraph/)

> a bioinformatic toolkit to align large sets of closely related genomes into a graph data structure

## Overview

**pangraph** provides both a command line interface, as well as a Julia library, to find homology amongst large collections of closely related genomes.
The core of the algorithm partitions each genome into _pancontigs_ that represent a sequence interval related by vertical descent.
Each genome is then an ordered walk along _pancontigs_; the collection of all genomes form a graph that captures all observed structural diversity.
**pangraph** is a standalone tool useful to parsimoniously infer horizontal gene transfer events within a community; perform comparative studies of genome gain, loss, and rearrangement dynamics; or simply to compress many related genomes.

## Installation

The core algorithm and command line tools are self-contained and require no additional dependencies.

### Library

Clone the repository
```bash
    git clone https://github.com/nnoll/pangraph.git && cd pangraph
```

Build the package. This will create a separate Julia environment for **pangraph**
```bash
    julia -e 'using Pkg; Pkg.build(".")'
```

Enter the REPL
```bash
    julia --project=.
```

### Relocatable binary
Releases can be obtained from [Github](https://github.com/nnoll/pangraph/releases)

## Examples

Align a multi-fasta `sequence.fa` and realign each _pancontig_ with MAFFT
```bash
	pangraph build sequence.fa | pangraph polish > graph.json
```

Export a graph `graph.json` to GFA for visualization
```bash
	pangraph export graph.json
```

Compute all pairwise graphs and estimate parsimonious number of events between strains.
Output all computed data to directory `pairs`
```bash
	pangraph marginalize -d pairs graph.json
```

See [Makefile](Makefile) for more examples.

## Citing
TBA

## License

[MIT License](LICENSE)
