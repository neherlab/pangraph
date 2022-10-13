# Updating the documentation

The documentation is hosted on <https://neherlab.github.io/pangraph/>.

### Automated builds

The Continuous integration (CI) will trigger a build and deployment of the docs website to GitHub Pages on every commit to `master` branch (a direct push or a merge of a pull request). You can track the build process on GitHub Actions:

https://github.com/neherlab/pangraph/actions

The docs CI is configured in

```bash
.github/workflows/docs.yml
```

### Build manually

To build and release a new version of the documentation manually:

1. run `make documentation`. This will create the directory `docs/build`
2. switch to the `gh-pages` branch of the repository
3. Substitute the content of the previously-created `docs/build` into the repo main directory and commit the changes.
