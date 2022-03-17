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
