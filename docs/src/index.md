# PanGraph
A **fast**, **self-contained** Julia library and command line tool suite to align multiple genomes into a pangenome graph.

## Introduction

Microbes don't evolve strictly by vertical descent and modification; in addition, there is persistent lateral exchange of genetic material within their local spatial neighborhood.
At the present, explanatory population genetic models exist at both extrema of limits horizontal transfer rate; researchers have a reasonable, quantitative understanding of both asexual and the fully decoupled, single-site, evolution.
However, there are currently no satisfactory models, and thus computational tools, in the more realistic intermediate regime, let alone that models which account for _any_ structural variation within a population.
In order to investigate such questions empirically in natural populations, we have developed a scalable multiple genome alignment tool, **PanGraph**, that identifies regions of mutual homology between large sets of closely related genomes.
This is expected to be useful to parsimoniously infer horizontal gene transfer events within a community; perform comparative studies of genome gain, loss, and rearrangement dynamics; or simply to compress many related genomes.

**PanGraph** (short for pangenome graph) builds a graph-based coordinate system for polymorphic variation in a microbial population that generalizes classical linear alignment techniques.
The graph is built by progressively aligning pairs of genomes using Minimap2's minimizer-based algorithm, and thus simultaneously aligns all genomes within the population in linear time.
The resultant graph represents contiguous intervals of homologous DNA as vertices and every genome as an ordered walk across such vertices.
Edges of the graph are unordered and only exist if at least one genome was found to connect both vertices in either the forward or reverse strand.
The documentation, and source code, uses the following terminology:

1. **Block/PanContig:**
    A contiguous interval of sequence that compresses a linear multiple sequence alignment and represents a single vertex within a pangraph.
    It is summarized by a consensus sequence, along with the polymorphisms needed for each genome to reconstruct the true sequence.
2. **Node:**
    An oriented traversal of a genome through a _block_, i.e. whether a given genome/duplication traverses a block along the forward or reverse direction.
    One genome can have multiple nodes that point to the same _block_.
3. **Path:**
    A sequence of _nodes_ that represent an ordered walk on the graph that reconstructs an input genome.
    Can be circular or linear.
4. **Edge:**
    A juxtaposition/breakpoint between two _blocks_ found in at least one genome.
5. **Graph/PanGraph:**
    A collection of _blocks_, associated to all recognized intervals of homology, and _paths_, genomes stored as an ordered walk of _nodes_.


## Installation

Pangraph is available:
- as a **Julia library**
- as a **Docker container**
- can be compiled into a **binary**


### as a Julia library

The library is written in and thus requires [Julia](https://julialang.org/downloads/) to be installed on your machine.

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

Note that to access the complete set of functionalities, the [optional dependencies](https://neherlab.github.io/pangraph/#Optional-dependencies) must be installed and available in your `$PATH`.


### Using Docker

Docker container image for PanGraph is available on Docker Hub: <https://hub.docker.com/r/neherlab/pangraph>

 **1. Install Docker**

Install Docker as described on the official website: <https://docs.docker.com/get-docker/>

Optionally setup Docker so that it runs without `sudo` on Linux: <https://docs.docker.com/engine/install/linux-postinstall/>

 **2. Pull a version of the PanGraph Docker image**

To obtain the latest released version, run:

```bash
    docker pull neherlab/pangraph:latest
```

To obtain a specific version, for example `1.2.3`, run:
   
```bash
    docker pull neherlab/pangraph:1.2.3
```

To obtain the latest development version (from `master` branch), run:

```bash
    docker pull neherlab/pangraph:master
```

> ⚠️ Note that the development versions can contain new and undocumented features, breaking changes and bugs. For most users, we recommend using either `:latest` or an explicit version.

**3. Run PanGraph Docker container**

Issue `docker run` command:

```bash
    docker run --rm -it \
      --name "pangraph-$(date +%s)" \
      --volume="$(pwd):/workdir" \
      --user="$(id -u):$(id -g)" \
      --workdir=/workdir \
      neherlab/pangraph:latest \
      bash -c "pangraph build --circular --alpha 0 --beta 0 /workdir/example_datasets/ecoli.fa.gz > graph.json"
```

Replace the `:latest` tag with either an explicit version, e.g. `:1.2.3` or `:master`, depending on which version you pulled in the previous section. If you haven't run `docker pull`, the `docker run` command should pull the corresponding version for you.

Here we mount current directory `.` (expressed as absolute path, using `pwd` shell command) as `/workdir` into the container so that pangraph can read the local
file `./example_datasets/ecoli.fa.gz` as `/workdir/example_datasets/ecoli.fa.gz"`:
    
```
                                 . -> /workdir
    ./example_datasets/ecoli.fa.gz -> /workdir/example_datasets/ecoli.fa.gz
```

The `--name` flag sets the name of the container and the `date` command there ensures that a unique name is created on every run. This is optional. The `--rm` flag deletes the container (but not the image) after run.

Replace `:latest` with a specific version if desired. The `:latest` tag can also be omitted, as it is the default. 


### Building binaries locally

PanGraph can be built locally on your machine by running (inside the cloned repo)
```bash
    export jc="path/to/julia/executable" make pangraph && make install
```
This will build the executable and place a symlink into `bin/`.

**Importantly,** if `jc` is not explicitly set, it will default to `vendor/julia-$VERSION/bin/julia`. If this file does not exist, we will download automatically for the user, provided the host system is Linux or MacOSX.
Moreover, for the compilation to work, it is necessary to have [MAFFT](https://mafft.cbrc.jp/alignment/software/) and [mmseqs2](https://github.com/soedinglab/MMseqs2) available in your `$PATH`, see [optional dependencies](https://neherlab.github.io/pangraph/#Optional-dependencies).

**Note,** it is [recommended by the PackageCompiler.jl documentation](https://julialang.github.io/PackageCompiler.jl/stable/#Installation-instructions) to utilize the officially distributed binaries for Julia, not those distributed by your Linux distribution. As such, compilation may not work if you attempt to do so.


### Optional dependencies

There are a few **optional** external programs that PanGraph can utilize[^1]:
1. [Mash](https://github.com/marbl/Mash) can be used to construct a guide tree in place of our internal algorithm (see [build](https://neherlab.github.io/pangraph/cli/build/) command options).
2. [MAFFT](https://mafft.cbrc.jp/alignment/software/) can be optionally used to polish block alignments (see [polish](https://neherlab.github.io/pangraph/cli/polish/) command). Only recommended for short alignments. 
3. [mmseqs2](https://github.com/soedinglab/MMseqs2) can be used as an alternative alignment kernel to the default *minimap2* (see [build](https://neherlab.github.io/pangraph/cli/build/) command options). It allows merging of more diverged sequences, at the cost of higher computational time.
4. [fasttree](http://www.microbesonline.org/fasttree/) is used to build phylogenetic trees for export in [PanX](https://github.com/neherlab/pan-genome-analysis)-compatible format (see [export](https://neherlab.github.io/pangraph/cli/export/) command options and the [tutorial section](https://neherlab.github.io/pangraph/tutorials/tutorial_3/#Explore-block-alignments-with-the-panX-visualization)).

In order to invoke all functionalities from PanGraph, these tools must be installed and available on `$PATH`.

If [conda](https://docs.conda.io/en/latest/) is available, one can run the following command to install all of these dependencies in a new environment named `pangraph`:
```
conda create -n pangraph -c conda-forge -c bioconda mmseqs2=13.45111 mash=2.2.2 mafft=7.475 FastTree=2.1.11
```

Alternatively, a script `bin/setup-pangraph` is provided within the repository to install both dependencies for a Linux machine without access to root. It assumes GNU coreutils are available.

These dependencies are already available within the Docker container.

[^1]: We recommend `mmseqs` version `13-45111`, `mash` version `v2.2.2`, `MAFFT` version `v7.475` and `fasttree` version `2.1.11`

## User's Guide

Basic functionality of **PanGraph** is provided by a command line interface.
This includes multiple genome alignment, the export of a genome alignment to various visualization formats, alignment polishing, and genome comparison tool. For more details please refer to the **Tutorials** section of the documentation.

Multithreading support is baked into the provided binary. Unfortunately, due to limitations in *julia*, the number of threads is set by the environment variable `JULIA_NUM_THREADS`

For uncovered use cases, functionality can be added by utilizing the underlying library functions.
Please see the high-level overview for definitions of library terminology.

## Citing PanGraph

PanGraph: scalable bacterial pan-genome graph construction. Nicholas Noll, Marco Molari, Liam P. Shaw, Richard Neher bioRxiv 2022.02.24.481757; doi: <https://doi.org/10.1101/2022.02.24.481757>
