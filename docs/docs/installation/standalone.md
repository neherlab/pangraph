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

|         | x86_64                                                                                                                                                                                                               | arm64                                                                                                                                                                                                                  |
|---------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| Linux   | [gnu](https://github.com/neherlab/pangraph/releases/latest/download/pangraph-x86_64-unknown-linux-gnu), [musl](https://github.com/neherlab/pangraph/releases/latest/download/pangraph-x86_64-unknown-linux-musl)[^1] | [gnu](https://github.com/neherlab/pangraph/releases/latest/download/pangraph-aarch64-unknown-linux-gnu), [musl](https://github.com/neherlab/pangraph/releases/latest/download/pangraph-aarch64-unknown-linux-musl)[^1] |
| macOS   | [download](https://github.com/neherlab/pangraph/releases/latest/download/pangraph-x86_64-apple-darwin)                                                                                                               | [download](https://github.com/neherlab/pangraph/releases/latest/download/pangraph-aarch64-apple-darwin)                                                                                                                |
| Windows | [download](https://github.com/neherlab/pangraph/releases/latest/download/pangraph-x86_64-pc-windows-gnu.exe)                                                                                                         | -                                                                                                                                                                                                                      |

:::info[Note]

There are 2 flavors of pangraph executables for Linux platforms - "gnu" and "musl". 

The "gnu" flavor is generally faster but requires GNU libc >= 2.17 to be present on the system. Linux distros like Debian, Ubuntu, Red Hat, Fedora, Arch all come with GNU libc by default. 

The "musl" flavor can be slower, but it makes no assumptions about libc. It is suitable for distros like Alpine, which have no GNU libc, or on older distros which have very old GNU libc.

All flavors require Linux kernel >= 3.2.

:::


:::info[Note]

On Unix platforms, don't forget to set executable file permission. You probably want to rename the executable and place it into `$PATH` to simplify its usage:

```bash
# Assuming ${HOME}/bin/pangraph is in $PATH
mv pangraph-x86_64-unknown-linux-gnu $HOME/bin/pangraph
chmod +x ${HOME}/bin/pangraph
pangraph --help
```

:::

## Installing from command line

You can also download pangraph programmatically, for example to add it to your Docker image or a bioinformatics pipeline:

<Tabs>
  <TabItem value="linuxX86" label="Linux x86">
    #### GNU, Latest Version
    ```bash
    curl -fsSL "https://github.com/neherlab/pangraph/releases/latest/download/pangraph-x86_64-unknown-linux-gnu" -o "pangraph" && chmod +x pangraph
    ```
    #### GNU, Specific Version
    ```bash
    curl -fsSL "https://github.com/neherlab/pangraph/releases/download/1.0.0/pangraph-x86_64-unknown-linux-gnu" -o "pangraph" && chmod +x pangraph
    ```
    #### MUSL, Latest Version
    ```bash
    curl -fsSL "https://github.com/neherlab/pangraph/releases/latest/download/pangraph-x86_64-unknown-linux-musl" -o "pangraph" && chmod +x pangraph
    ```
    #### MUSL, Specific Version
    ```bash
    curl -fsSL "https://github.com/neherlab/pangraph/releases/download/1.0.0/pangraph-x86_64-unknown-linux-musl" -o "pangraph" && chmod +x pangraph
    ```
  </TabItem>
  <TabItem value="linuxArm64" label="Linux ARM64">
    #### GNU, Latest Version
    ```bash
    curl -fsSL "https://github.com/neherlab/pangraph/releases/latest/download/pangraph-aarch64-unknown-linux-gnu" -o "pangraph" && chmod +x pangraph
    ```
    #### GNU, Specific Version
    ```bash
    curl -fsSL "https://github.com/neherlab/pangraph/releases/download/1.0.0/pangraph-aarch64-unknown-linux-gnu" -o "pangraph" && chmod +x pangraph
    ```
    #### MUSL, Latest Version
    ```bash
    curl -fsSL "https://github.com/neherlab/pangraph/releases/latest/download/pangraph-aarch64-unknown-linux-musl" -o "pangraph" && chmod +x pangraph
    ```
    #### MUSL, Specific Version
    ```bash
    curl -fsSL "https://github.com/neherlab/pangraph/releases/download/1.0.0/pangraph-aarch64-unknown-linux-musl" -o "pangraph" && chmod +x pangraph
    ```
  </TabItem>
  <TabItem value="macOSx86" label="macOS x86">
    #### Latest Version
    ```bash
    curl -fsSL "https://github.com/neherlab/pangraph/releases/latest/download/pangraph-x86_64-apple-darwin" -o "pangraph" && chmod +x pangraph
    ```
    #### Specific Version
    ```bash
    curl -fsSL "https://github.com/neherlab/pangraph/releases/download/1.0.0/pangraph-x86_64-apple-darwin" -o "pangraph" && chmod +x pangraph
    ```
  </TabItem>
  <TabItem value="macOSArm64" label="macOS ARM64">
    #### Latest Version
    ```bash
    curl -fsSL "https://github.com/neherlab/pangraph/releases/latest/download/pangraph-aarch64-apple-darwin" -o "pangraph" && chmod +x pangraph
    ```
    #### Specific Version
    ```bash
    curl -fsSL "https://github.com/neherlab/pangraph/releases/download/1.0.0/pangraph-aarch64-apple-darwin" -o "pangraph" && chmod +x pangraph
    ```
  </TabItem>
  <TabItem value="windowsX86" label="Windows x86">
    #### Latest Version
    ```bash
    curl -fsSL "https://github.com/neherlab/pangraph/releases/latest/download/pangraph-x86_64-pc-windows-gnu.exe" -o "pangraph.exe" && chmod +x pangraph.exe
    ```
    #### Specific Version
    ```bash
    curl -fsSL "https://github.com/neherlab/pangraph/releases/download/1.0.0/pangraph-x86_64-pc-windows-gnu.exe" -o "pangraph.exe" && chmod +x pangraph.exe
    ```
  </TabItem>
</Tabs>





## optional dependencies - mmseqs2

When building graphs, in addition to the native _minimap2_, Pangraph can use [_mmseqs2_](https://github.com/soedinglab/MMseqs2) as an alignment kernel.[^1] _mmseqs2_ is able to find homology between sequences with up to around 30% sequence divergence, at the cost of higher computational time.

Use of this kernel requires installation of _mmseqs2_. This can be done for example using conda with:

```bash
conda install -c conda-forge -c bioconda mmseqs2
```

:::note

Pangraph was tested with _mmseqs2_ version 17.b804f.

:::

[^1]: This kernel can be selected with the `-k mmseqs` option when running the [`pangraph build`](../reference.md#pangraph-build) command.
