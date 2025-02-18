# PanGraph

[![Documentation](https://img.shields.io/badge/Documentation-Link-blue.svg)](https://docs.pangraph.org/)
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

The original implementation of pangraph (version v0) was implemented in Julia and was described in the publication [Noll, Molari, Shaw and Neher, 2023](https://www.microbiologyresearch.org/content/journal/mgen/10.1099/mgen.0.001034).
The current version (v1) is a reimplementation of the original algorithm in Rust by Ivan Aksamentov and Marco Molari.
The new implementation should be much easier to install and is faster in many use cases.

## Installation

Pangraph is available:
- as a **standalone binary**
- as a **docker container**

For more extended instructions on installation please refer to the [documentation](https://docs.pangraph.org/category/installation).

### Standalone binary

This is the recommended way to install Pangraph. You can download the latest release for your operating system [from here](https://docs.pangraph.org/installation/standalone).

### Docker container

PanGraph is available as a Docker container:

```bash
docker pull neherlab/pangraph:latest
```

See the [documentation](https://docs.pangraph.org/installation/with-docker) for extended instructions on its usage.


## Examples

Please refer to the [tutorials within the documentation](https://docs.pangraph.org/category/tutorial) for an in-depth usage guide.
For a quick reference, see below.

Align a multi-fasta `sequences.fa` in a graph:
```bash
pangraph build sequences.fa -o graph.json
```

Extract the core-genome alignment from the graph, with blocks appearing in the order of the reference genome `NC_010468`:
```bash
pangraph export core-genome graph.json \
  --guide-strain NC_010468 \
  -o core_genome_aln.fa
```

Export the graph in gfa format for visualization:
```bash
pangraph export gfa graph.json -o graph.gfa
```

Reconstruct input sequences from the graph:
```bash
pangraph reconstruct graph.json -o sequences.fa
```

## PyPangraph

PyPangraph is a python package with convenient utilities to load and explore the graph data structure, see the [documentation](https://docs.pangraph.org/pypangraph/about-pypangraph) for installation instructions and more examples.

```python
import pypangraph as pp

graph = pp.Pangraph.load_graph("graph.json")
print(graph)
# pangraph object with 15 paths, 137 blocks and 1042 nodes
```


## Citation
If you use PanGraph in scientific publications, please cite the original paper presenting the algorithm:

*PanGraph: scalable bacterial pan-genome graph construction*
Nicholas Noll, Marco Molari, Liam P. Shaw, Richard A. Neher
Microbial Genomics, **9**(6), 2023; doi: [10.1099/mgen.0.001034](https://doi.org/10.1099/mgen.0.001034)


## License

[MIT License](LICENSE)

> [!NOTE]
> The legacy v0 version of Pangraph is now stored on the [`v0` branch](https://github.com/neherlab/pangraph/tree/v0) of the repository, and legacy documentation is available [here](https://v0.docs.pangraph.org/).
