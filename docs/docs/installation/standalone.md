---
sidebar_position: 1
---

import Tabs from '@theme/Tabs';
import TabItem from '@theme/TabItem';

# Installation: standalone

:::info[Note]

This is the recommended way to install Pangraph

:::

Pangraph CLI is a self-contained single-file executable which can be downloaded and run without any particular setup.



|         | x86_64                                                                                                                                                                                                            | arm64                                                                                                                                                                                                               |
| ------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Linux   | [gnu](https://github.com/neherlab/pangraph/releases/latest/download/pangraph-x86_64-unknown-linux-gnu), [musl](https://github.com/neherlab/pangraph/releases/latest/download/pangraph-x86_64-unknown-linux-musl)* | [gnu](https://github.com/neherlab/pangraph/releases/latest/download/pangraph-aarch64-unknown-linux-gnu), [musl](https://github.com/neherlab/pangraph/releases/latest/download/pangraph-aarch64-unknown-linux-musl)* |
| macOS   | [download](https://github.com/neherlab/pangraph/releases/latest/download/pangraph-x86_64-apple-darwin)                                                                                                            | [download](https://github.com/neherlab/pangraph/releases/latest/download/pangraph-aarch64-apple-darwin)                                                                                                             |
| Windows | [download](https://github.com/neherlab/pangraph/releases/latest/download/pangraph-x86_64-pc-windows-gnu.exe)                                                                                                      | -                                                                                                                                                                                                                   |
## optional dependencies - mmseqs2

When building graphs, in addition to the native _minimap2_, Pangraph can use [_mmseqs2_](https://github.com/soedinglab/MMseqs2) as an alignment kernel.[^1] _mmseqs2_ is able to find homology between sequences with up to around 30% sequence divergence, at the cost of higher computational time.

[^1]: This kernel can be selected with the `-k mmseqs` option when running the [`pangraph build`](../reference.md#pangraph-build) command.

Use of this kernel requires installation of _mmseqs2_. This can be done for example using conda with:

```bash
conda install -c conda-forge -c bioconda mmseqs2
```

:::note

Pangraph was tested with _mmseqs2_ version 17.b804f.

:::