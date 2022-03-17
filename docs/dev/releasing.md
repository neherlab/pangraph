## 🆕 Releasing pangraph

### Releasing a new version

Continuous integration (CI) will build a new version of the Docker container (see `Dockerfile`) on every pushed git tag.

In order to creat and push a git tag, run:

```
git tag $RELEASE_VERSION
git push origin --tags
```

where `$RELEASE_VERSION` is a valid [semantic version](https://semver.org/), without a `v` prefix (i.e. `1.2.3` is correct, `v1.2.3` is not).


### Monitoring and debugging CI build

The status of the builds can be seen on GitHub Actions page:

https://github.com/neherlab/pangraph/actions

### Verifying CI build

After CI build successfully finishes, pull and run the new version to make sure it works as expected


```bash
docker pull neherlab/pangraph:$RELEASE_VERSION

docker run --rm -it \
  --name "pangraph-$(date +%s)" \
  --volume="$(pwd)/path-to-fasta:/workdir" \
  --workdir=/workdir pangraph:$RELEASE_VERSION \
  bash -c "pangraph build --circular --alpha 0 --beta 0 /workdir/test.fa"
```

Here we mount local directory `path-to-fasta` as `/workdir` so that pangraph can read the `/workdir/test.fa"` file. 

> 👷 TODO: implement automated tests


### Modifying continuous integration workflow

See `.github/workflows/build.yml`


### Modifying Docker image

See `Dockerfile`
