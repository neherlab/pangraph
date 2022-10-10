# pangraph

[![Documentation](https://img.shields.io/badge/Documentation-Link-blue.svg)](https://neherlab.github.io/pangraph/)
![Docker Image Version (latest semver)](https://img.shields.io/docker/v/neherlab/pangraph?label=docker)
![Docker Pulls](https://img.shields.io/docker/pulls/neherlab/pangraph)

> a bioinformatic toolkit to align large sets of closely related genomes into a graph data structure

## Overview

**pangraph** provides both a command line interface, as well as a Julia library, to find homology amongst large collections of closely related genomes.
The core of the algorithm partitions each genome into _pancontigs_ that represent a sequence interval related by vertical descent.
Each genome is then an ordered walk along _pancontigs_; the collection of all genomes form a graph that captures all observed structural diversity.
**pangraph** is a standalone tool useful to parsimoniously infer horizontal gene transfer events within a community; perform comparative studies of genome gain, loss, and rearrangement dynamics; or simply to compress many related genomes.

## Installation

The core algorithm and command line tools are self-contained and require no additional dependencies.
The library is written in and thus requires [Julia](https://julialang.org/downloads/) to be installed on your machine.

### Julia Library

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
    julia -e 'using Pkg; Pkg.add(url="https://github.com/nnoll/minimap2_jll.jl"); Pkg.add(url="https://github.com/neherlab/pangraph.git")'
```

The PanGraph package is available globally within the Julia REPL.

### Relocatable binary

**pangraph** can be built locally on your machine by running (inside the cloned repo)
```bash
    export jc="path/to/julia/executable" make pangraph && make install
```
This will build the executable and place a symlink into `bin/`.
**Importantly,** if `jc` is not explicitly set, it will default to `vendor/julia-$VERSION/bin/julia`. If this file does not exist, we will download automatically for the user, provided the host system is Linux or MacOSX.
Moreover, for the compilation to work, it is necessary to have [MAFFT](https://mafft.cbrc.jp/alignment/software/) and [mmseqs2](https://github.com/soedinglab/MMseqs2) available in your `$PATH`, see [optional dependencies](#optional-dependencies).
**Note,** it is [recommended by the PackageCompiler.jl documentation](https://julialang.github.io/PackageCompiler.jl/stable/#Installation-instructions) to utilize the officially distributed binaries for Julia, not those distributed by your Linux distribution. As such, compilation may not work if you attempt to do so.

### Optional dependencies
**pangraph** can _optionally_ use [mash](https://github.com/marbl/Mash), [MAFFT](https://mafft.cbrc.jp/alignment/software/) or [mmseqs2](https://github.com/soedinglab/MMseqs2).
For full functionality, it is recommended to install these tools and have them available on `$PATH`.

Alternatively, a script `bin/setup-pangraph` is provided to install both tools into `bin/` for Linux-based operating systems.

### Docker container

PanGraph is available as a Docker container:

```bash
    docker pull neherlab/pangraph:latest
```

See the [documentation](https://neherlab.github.io/pangraph/#Installation) for more extended instuctions on its usage.

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

## Citing
PanGraph: scalable bacterial pan-genome graph construction
Nicholas Noll, Marco Molari, Richard Neher
bioRxiv 2022.02.24.481757; doi: https://doi.org/10.1101/2022.02.24.481757

## License

[MIT License](LICENSE)
