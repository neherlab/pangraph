# Updating the documentation

The documentation is hosted on <https://neherlab.github.io/pangraph/>.

### Automated documentation releases

The Continuous integration (CI) will trigger a build and deployment of the docs website to GitHub Pages on every commit to `master` branch (a direct push or a merge of a pull request). You can track the build process on GitHub Actions:

https://github.com/neherlab/pangraph/actions

The docs CI is configured in

```bash
.github/workflows/docs.yml
```

### Build manually

To build a new version of the documentation manually, run `make documentation`. This will create the directory `docs/build`. 

### Preview locally

You can preview the resulting docs website locally by serving the resulting `docs/build`, using any static file server. For example, if you have [Node.js](https://nodejs.org/en/download/) and npx installed, you could run `serve` package:

```bash
npx serve --listen=tcp://localhost:8888 docs/build
```

and then open `http://localhost:8888` in a browser.

### Release manually

To release a new version of the documentation manually (bypassing the CI workflow):

1. build the documentation 
2. switch to the `gh-pages` branch of the repository
3. Substitute the content of the previously-created `docs/build` into the repo main directory and commit the changes.
