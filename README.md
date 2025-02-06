# pangraph

[![Documentation](https://img.shields.io/badge/Documentation-Link-blue.svg)](https://neherlab.github.io/pangraph/)
![Docker Image Version (latest semver)](https://img.shields.io/docker/v/neherlab/pangraph?label=docker)
![Docker Pulls](https://img.shields.io/docker/pulls/neherlab/pangraph)

> a bioinformatic toolkit to align large sets of closely related genomes into a graph data structure

> [!WARNING]  
> Pangraph is currently undergoing a major migration between v0 and v1. In this short transition period links and documentation may be inconsistent.

## Overview

**pangraph** provides a command line interface to find homology amongst large collections of closely related genomes.
The core of the algorithm partitions each genome into _blocks_ that represent a sequence interval related by vertical descent.
Each genome is then an ordered walk along _blocks_. The collection of all genomes form a graph that captures all observed structural diversity.
**pangraph** is a standalone tool useful to parsimoniously infer horizontal gene transfer events within a community; perform comparative studies of genome gain, loss, and rearrangement dynamics; or simply to compress many related genomes.


## Installation

Pangraph is available:
- as a **standalone binary**
- as a **conda package**
- as a **docker container**

For more extended instructions on installation please refer to the documentation.

### Standalone binary

### Conda package

### Docker container

PanGraph is available as a Docker container:

```bash
    docker pull neherlab/pangraph:latest
```

See the documentation for extended instuctions on its usage.


## Examples

Please refer to the tutorials within the documentation for an in-depth usage guide.
For a quick reference, see below.

Align a multi-fasta `sequences.fa` in a graph:
```bash
	pangraph build sequences.fa > graph.json
```

Extract the core-genome alignment from the graph, with blocks appearing in the order of the reference genome `NC_010468`:
```bash
	pangraph export core-genome graph.json \
        --guide-strain NC_010468 \
        > core_genome_aln.fa
```

Export the graph in gfa format for visualization:
```bash
    pangraph export gfa graph.json > graph.gfa
```

Reconstruct input sequences from the graph:
```bash
    pangraph reconstruct graph.json > sequences.fa
```

## PyPangraph

PyPangraph is a python package with convenient utilities to load and explore the graph data structure, see the documentation for more details.

```python
import pypangraph as pp

graph = pp.Pangraph.load_graph("graph.json")
print(graph)
# pangraph object with 15 paths, 137 blocks and 1042 nodes
```


## Citing
PanGraph: scalable bacterial pan-genome graph construction
Nicholas Noll, Marco Molari, Richard Neher
bioRxiv 2022.02.24.481757; doi: https://doi.org/10.1101/2022.02.24.481757


## License

[MIT License](LICENSE)
