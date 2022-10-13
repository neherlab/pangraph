# Releasing

### Releasing a new version

Continuous integration (CI) will build a new version of the Docker container (see `Dockerfile`) on every pushed git tag.

Make sure you are on a correct branch and commit. Most of the time you want to release code from `master`:

```bash
git checkout master
```

In order to create and push a git tag, run:

```
git tag $RELEASE_VERSION
git push origin --tags
```

where `$RELEASE_VERSION` is a valid [semantic version](https://semver.org/), without a `v` prefix (i.e. `1.2.3` is correct, `v1.2.3` is not).

The CI workflow will build the container image and will push it to Docker Hub. The image will be tagged with:
 
 - `latest` (and will overwrite existing `latest` tag there)
 - `$RELEASE_VERSION`

Both tags should point to the same image, i.e. their sha hashes should be exactly the same.

This image version can then be referred to as:

 - `neherlab/pangraph:$RELEASE_VERSION`
 - `neherlab/pangraph:latest`
 - `neherlab/pangraph` (which is the same as `neherlab/pangraph:latest`)

for example in `docker pull` and `docker run` commands.


### Monitoring and debugging CI build

The status of the builds can be seen on GitHub Actions page:

https://github.com/neherlab/pangraph/actions

### Verifying CI build

After CI build successfully finishes, check Docker Hub to ensure that the new tag is present and that the "latest" tag is updated and points to the same hash:

https://hub.docker.com/r/neherlab/pangraph

Pull and run the new version to make sure it works as expected:

```bash
docker pull neherlab/pangraph:$RELEASE_VERSION

docker run --rm -it \
  --name "pangraph-$(date +%s)" \
  --volume="$(pwd)/path-to-fasta:/workdir" \
  --user="$(id -u):$(id -g)" \
  --workdir=/workdir neherlab/pangraph:$RELEASE_VERSION \
  bash -c "pangraph build --circular --alpha 0 --beta 0 /workdir/test.fa"
```

Here we mount local directory `path-to-fasta` as `/workdir` so that pangraph can read the `/workdir/test.fa"` file. 

> 👷 TODO: implement automated tests


### Modifying continuous integration workflow

See `.github/workflows/build.yml`


### Modifying Docker image

See `Dockerfile`
