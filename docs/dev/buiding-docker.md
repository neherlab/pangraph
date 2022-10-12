## 👷 Building pangraph Docker image locally

### Install dependencies

 - Install bash and make.

 - Install Docker: https://docs.docker.com/get-docker/

 - Optionally setup Docker so that it runs without `sudo`: https://docs.docker.com/engine/install/linux-postinstall/


### Build Docker image locally

Run:

```bash
make docker
```

This will build the Docker image tagged `neherlab/pangraph` (more precisely `neherlab/pangraph:latest`). If already exists, it will replace the existing image with that tag. The build will take some time.

If completed successfully, then the image can be used right away. Refer to user documentation. Skip the "Pull Docker image" step.


### Explore contents, layers and optimize image size

You could use [dive tool](https://github.com/wagoodman/dive) to see what's inside an image:

```bash
dive neherlab/pangraph:<tag>
```

Each [layer](https://stackoverflow.com/questions/31222377/what-are-docker-image-layers) reflects `FROM`, `COPY` and `RUN` commands and the files that have been added to the overlay file system of the image. This can be used to find redundant files. You could then further optimize `Dockerfile` and make the image smaller.


### Testing the image

Once the image has been built locally it can be tested by running (from the repository main folder):

```bash
docker run -it --rm \
    --volume="$(pwd):/workdir" \
    --workdir="/workdir" \
    --user="$(id -u):$(id -g)" \
    neherlab/pangraph \
    bash docs/dev/docker_test.sh
```

Or more simply using instructions in the Makefile:

```bash
make docker-test
```

This will test all the available commands, see `docs/dev/docker_tests.sh` script.
