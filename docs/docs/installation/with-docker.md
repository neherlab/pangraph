---
sidebar_position: 2
---

# With docker

Run Pangraph CLI inside a Docker container.

:::info[Note]

The following instructions assume you have [Docker](https://www.docker.com/) installed

:::

Pull latest version of [`neherlab/pangraph` image from DockerHub](https://hub.docker.com/r/neherlab/pangraph):

```shell
docker pull neherlab/pangraph
```

or pull a particular version

```shell
docker pull neherlab/pangraph:1.0
```

Run pangraph commands described in [Usage](../reference.md) as follows:

```shell
docker run -it --rm neherlab/pangraph --help
docker run -it --rm neherlab/pangraph build my.fasta
```


