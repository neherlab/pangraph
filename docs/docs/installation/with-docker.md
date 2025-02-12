---
sidebar_position: 2
---

# Installation: with docker

Run Pangraph CLI inside a Docker container.

:::info[Note]

The following instructions assume you have [Docker](https://www.docker.com/) installed and is familiar with the basics. Explanation of Docker is out of scope of this documentation. Refer to [Docker documentation](https://docs.docker.com) if necessary.

:::

Pull latest version of [`neherlab/pangraph` image from DockerHub](https://hub.docker.com/r/neherlab/pangraph):

```shell
docker pull neherlab/pangraph
```

or pull a particular version

```shell
docker pull neherlab/pangraph:1
```

Run pangraph commands described in [CLI reference](../reference.md) or [Tutorials](/category/tutorial) as in following examples:

```shell
docker run -it --rm neherlab/pangraph pangraph --help

# Mount current directory to /data and run as current user
docker run -it --rm --volume="$(pwd):/data" --user="$(id -u):$(id -g)" neherlab/pangraph pangraph build /data/input.fasta -o /data/pangraph.json
```

:::info[Note]

Don't forget to [add volumes](https://docs.docker.com/engine/storage/volumes/#syntax) to be able to access data on the host. You might also want to [run docker container as a non-root user](https://docs.docker.com/engine/containers/run/#user).

:::


## Extending the image

You can base your own containers based on pangraph image and add your own software and functionality, for example:

```dockerfile
FROM neherlab/pangraph:1

RUN sudo apt-get update && sudo apt-get install python3

COPY my-script.py /usr/bin/my-script.py

# ... remainder of your image

```

or you can add [standalone](standalone.md) pangraph binaries to any image:

```dockerfile
FROM my-org/my-image

RUN curl -fsSL "https://github.com/neherlab/pangraph/releases/latest/download/pangraph-x86_64-unknown-linux-gnu" -o "/usr/bin/pangraph" && chmod +x /usr/bin/pangraph

# ... remainder of your image

```
