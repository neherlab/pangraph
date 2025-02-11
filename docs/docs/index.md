---
sidebar_position: 1
---

# Introduction to Pangraph

:::danger[Under construction 👷]

Pangraph is currently under heavy development. Bugs and crashes are to be expected.

:::


The content and structure of bacterial genomes evolves very rapidly:
Part of the genome can be cut out, duplicated, or inverted.
In addition, genomic material can be gained from the outside for example through phage infection or DNA uptake and integration.
As a result, comparing bacterial genomes is more complicated than analyzing differences in the alignment of homologous sequences.
Instead, one would like to understand how diversity in terms of content and structure has arisen through insertions, deletions, transpositions over the course of evolution.
To address such questions, we have developed a scalable multiple genome alignment tool, **PanGraph**, that identifies regions of mutual homology between large sets of closely related genomes and represents them in a graph.

![img](./assets/t1_main_scheme.png)

This is expected to be useful to parsimoniously infer horizontal gene transfer events within a community; perform comparative studies of genome gain, loss, and rearrangement dynamics; or simply to compress many related genomes.

The resultant graph represents contiguous intervals of homologous DNA as vertices and every genome as an ordered walk across such vertices.
Edges of the graph are unordered and only exist if at least one genome was found to connect both vertices in either the forward or reverse strand.
For a more detailed description of the graph structure, see [what is a pangraph](tutorial/t01-building-pangraph.md#what-is-a-pangraph).

This documentation is structures as a [set of tutorials](/category/tutorial) that explain the essential steps to build and manipulate a graph, along with a [reference documentation](/reference) of the available commands. In addition, we provide a python library [pyPanGraph](/category/pypangraph) for analysis of the graph data structure in Python.


