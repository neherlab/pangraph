# PanGraph
A **fast**, **self-contained** Julia library and command line tool suite that aligns multiple genomes simultaneously.

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

There are multiple ways to install PanGraph (either the library or just command line interface)

### From Julia REPL
```julia
    (@v1.x) pkg> add https://github.com/nnoll/pangraph.git
```

### From Command Line
```bash
    julia -e 'using Pkg; Pkg.add("https://github.com/nnoll/pangraph.git"); Pkg.build()'
```

### Local Environment

Clone the repository.
```bash
    git clone https://github.com/nnoll/pangraph.git && cd pangraph
```

Build the package. This will create a seperate Julia environment for PanGraph
```bash
    julia -e 'using Pkg; Pkg.build(".")'
```

Enter the REPL
```bash
    julia --project=.
```

### Binary
Additionally, `pangraph` is available as a standalone, relocatable binary that should work on any Linux or MacOSX machine.
Releases can be obtained from [Github](https://github.com/nnoll/pangraph/releases)

### Optional dependencies

There are a few **optional** external programs that PanGraph can utilize
1. [Mash](https://github.com/marbl/Mash) can be used to construct a guide tree in place of our internal algorithm.
2. [MAFFT](https://mafft.cbrc.jp/alignment/software/) can be optionally used to polish homologous alignments. Only recommended for short alignments.
In order to invoke functionality from PanGraph, these tools must be installed and available on **$PATH**

For convenience, a script `bin/setup-pangraph` is provided within the repository to install both dependencies for a Linux machine without access to root.
It assumes GNU coreutils are available.

## User's Guide

Basic functionality of PanGraph is provided by a command line interface.
This includes multiple genome alignment, the export of a genome alignment to various visualization formats, alignment polishing, and genome comparison tool.
Additionally, generation of basic synthetic data is included for testing.

For uncovered use cases, functionality can be added by utilizing the underlying library functions.
Please see the high-level overview for definitions of library terminology.

## Citing PanGraph

Add citation here