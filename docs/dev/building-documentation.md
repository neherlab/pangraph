# Updating the documentation

The documentation is hosted on <https://neherlab.github.io/pangraph/>.

To build and release a new version of the documentation:

1. run `make documentation`. This will create the directory `docs/build`
2. switch to the `gh-pages` branch of the repository
3. Substitute the content of the previously-created `docs/build` into the repo main directory and commit the changes.