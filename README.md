# pangraph

[![Documentation](https://img.shields.io/badge/Documentation-Link-blue.svg)](https://neherlab.github.io/pangraph/)

> a bioinformatic toolkit to align large sets of closely related genomes into a graph data structure

## Overview

**pangraph** provides both a command line interface, as well as a Julia library, to find homology amongst large collections of closely related genomes.
The core of the algorithm partitions each genome into _pancontigs_ that represent a sequence interval related by vertical descent.
Each genome is then an ordered walk along _pancontigs_; the collection of all genomes form a graph that captures all observed structural diversity.
**pangraph** is a standalone tool useful to parsimoniously infer horizontal gene transfer events within a community; perform comparative studies of genome gain, loss, and rearrangement dynamics; or simply to compress many related genomes.

## Installation

The core algorithm and command line tools are self-contained and require no additional dependencies.
The library is written in and thus requires Julia to be installed on your machine.
Julia binaries for all operating systems can be found [here](https://julialang.org/downloads/).

### Library

#### Local Environment

Clone the repository
```bash
    git clone https://github.com/neherlab/pangraph.git && cd pangraph
```

Build the package. This will create a separate Julia environment for **pangraph**
```bash
    julia --project=. -e 'using Pkg; Pkg.build()'
```

Enter the REPL
```bash
    julia --project=.
```

#### Global Package

**Important** please do not mix this method with that described above.
Instead of creating a _local_ PanGraph specific environment, this method will install into the Julia base environment.
We recommend, unless for a specific reason, to default to installing within a local environment.
However, if needed, global installation can be achieved by running

```bash
    julia -e 'using Pkg; Pkg.add(url="https://github.com:nnoll/minimap2_jll.jl"); Pkg.add(url="https://github.com:neherlab/pangraph.git")'
```

The PanGraph package is available globally within the Julia REPL.

### Relocatable binary
Releases can be obtained from [GitHub](https://github.com/neherlab/pangraph/releases)

Alternatively, **pangraph** can be built locally on your machine by running (inside the cloned repo)
```bash
    export jc="path/to/julia/executable" make pangraph && make install
```
This will build the executable and place a symlink into `bin/`.
**Importantly,** if `jc ` is not explicitly set, it will default to vendor/julia-$VERSION/bin/julia.
If this file does not exist, we will download automatically for the user, provided the host system is Linux or MacOSX.
**Note,** it is recommended by the PackageCompiler.jl documentation to utilize the officially distributed binaries, not those distributed by your Linux distribution.
As such, it may not work if you attempt to do so.

### Optional dependencies
**pangraph** can _optionally_ use both [mash](https://github.com/marbl/Mash) and [MAFFT](https://mafft.cbrc.jp/alignment/software/).
For full functionality, it is recommended to install these tools and have them available on `$PATH`.

Alternatively, a script `bin/setup-pangraph` is provided to install both tools into `bin/` for Linux-based operating systems.

## Examples

Please refer to the tutorials within the [documentation](https://neherlab.github.io/pangraph/) for an in-depth usage guide.
For a quick reference, see below.

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

See [Makefile](Makefile) for more real-world examples.

## Citing
TBA

## License

[MIT License](LICENSE)
