---
sidebar_position: 4
---

# Optional Dependencies

Pangraph can be used with additional tools to improve its functionality.
Below are some optional dependencies that can be installed to enable additional features.

## mmseqs2

When building graphs, in addition to the native _minimap2_, Pangraph can use [_mmseqs2_](https://github.com/soedinglab/MMseqs2) as an alignment kernel. _mmseqs2_ is able to find homology between sequences with up to around 30% sequence divergence, at the cost of higher computational time.

This kernel can be selected with the `-k mmseqs` option when running the [`pangraph build`](../usage/reference#pangraph-build) command.

To install _mmseqs2_ using conda, run the following command:

```bash
conda install -c conda-forge -c bioconda mmseqs2
```