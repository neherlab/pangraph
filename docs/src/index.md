# PanGraph
A **fast**, **self-contained** Julia library and command line tool suite that aligns multiple genomes simultaneously.

## Introduction

Add stuff here

## Installation

There are multiple ways to install PanGraph

### From Julia REPL
```julia
    (@v1.x) pkg> add https://github.com/nnoll/pangraph.git
```

### From Command Line
```bash
    julia -e 'using Pkg; Pkg.add("https://github.com/nnoll/pangraph.git")'
```

### Local Environment

Clone the repository.
```bash
    git clone https://github.com/nnoll/pangraph.git && cd pangraph
```

Instantiate the package. This will create a seperate Julia environment for PanGraph
```bash
    julia -e 'using Pkg; Pkg.instantiate(".")'
```

Enter the REPL
```bash
    julia --project=.
```

### Optional dependencies

There are a few **optional** external programs that PanGraph can utilize
1. [Mash](https://github.com/marbl/Mash) can be used to construct a guide tree in place of our internal algorithm.
2. [MAFFT](https://mafft.cbrc.jp/alignment/software/) can be optionally used to polish homologous alignments. Only recommended for short alignments.
In order to invoke functionality from PanGraph, these tools must be installed and available on **$PATH**

## User's Guide

## Citing PanGraph

Add citation here
