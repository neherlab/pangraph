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

pangraph is available:
- as a **julia library**
- as a **Docker container**
- it can be compiled into a relocatable **binary**

For more extended instructions on installation please refer to the [documentation](https://neherlab.github.io/pangraph/#Installation).

### Julia Library

To install pangraph as a julia library in a local environment:
```bash
    # clone the repository
    git clone https://github.com/neherlab/pangraph.git && cd pangraph
    # build the package
    julia --project=. -e 'using Pkg; Pkg.build()'
```

The library can be accessed directly by entering the REPL:
```bash
    julia --project=.
```

Alternatively, command-line functionalities can be accessed by running the main `src/PanGraph.jl` script:
```bash
    # example: build a graph from E.coli genomes
    julia --project=. src/PanGraph.jl build -c example_datasets/ecoli.fa.gz > graph.json
```

Note that to access the complete set of functionalities, the [optional dependencies](#optional-dependencies) must be installed and available in your `$PATH`.


### Docker container

PanGraph is available as a Docker container:

```bash
    docker pull neherlab/pangraph:latest
```

See the [documentation](https://neherlab.github.io/pangraph/#Installation) for extended instuctions on its usage.


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

**pangraph** can _optionally_ use [mash](https://github.com/marbl/Mash), [MAFFT](https://mafft.cbrc.jp/alignment/software/), [mmseqs2](https://github.com/soedinglab/MMseqs2) or [fasttree](http://www.microbesonline.org/fasttree/) for some optional functionalities, as explained in [the documentation](https://neherlab.github.io/pangraph/#Optional-dependencies).
For use of these functionalities, it is recommended to install these tools and have them available on `$PATH`.

Alternatively, a script `bin/setup-pangraph` is provided to install both tools into `bin/` for Linux-based operating systems.


## Examples

Please refer to the tutorials within the [documentation](https://neherlab.github.io/pangraph/) for an in-depth usage guide.
For a quick reference, see below.

Align a multi-fasta `sequence.fa` and realign each _pancontig_ with MAFFT
```bash
	pangraph build sequence.fa | pangraph polish > graph.json
```

Export a graph `graph.json` into `export/pangraph.gfa` as GFA for visualization
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
