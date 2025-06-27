# Documentation

## Writing docs

Documentation source code is under [`docs/`](../docs). It is built using [Docusaurus](https://docusaurus.io/).

### Requirements:

- Node.js. Download it from the official website: [nodejs.org](https://nodejs.org), or use [nvm](https://github.com/nvm-sh/nvm). The only supported version is the one that's declared in `.nvmrc` file. Different versions as well as installations from Linux package repositories, Homebrew, Conda etc. are not supported - they might work, but you are on your own.

- bun (https://bun.sh)

### Steps

This will install nvm, Node.js, bun, and will run development server of the docs website.

```bash
# Install nvm. Read this on how to add it to $PATH: https://github.com/nvm-sh/nvm
curl -fssLo- https://raw.githubusercontent.com/nvm-sh/nvm/master/install.sh | bash

# Install and use Node.js version from .nvmrc file
nvm install && nvm use

# Install bun
npm -g install bun

# Install npm dependencies
bun install

# Run dev server
bun dev
```

The development server will serve the documentation website on local port `4000`. Open [localhost:4000](http://localhost:4000) in a browser to see it. Edits to the content and code changes will be reflected on the fly.

The documentation files are written in markdown (`.md`) or (`.mdx`) format and are located under [`docs/docs/`](../docs/docs).

## Release build

If you want to build and inspect the production build of the docs, run:

```bash
bun prod:build
bun prod:serve
```

The `build` command will produce the directory with HTML, CSS and JS files ready to be deployed to a web hosting. The `serve` command will serve these files on port `5000`.

## Generate command-line reference documentation

The `docs/docs/reference.md` file is generated using script `generate-reference-docs`:

```bash
cargo build --bin=pangraph
cd docs
./generate-reference-docs "../target/debug/pangraph" "docs/reference.md"
```

> ⚠️ Do not edit the generated file manually! All manual changes will be overwritten by automation.

## Deployment

The documentation is automatically deployed on changes to `release-docs` branch. See [.github/workflows/docs.yml](../.github/workflows/docs.yml).

If you have documentation changes on `master` branch, which you want to release to the publicly facing documentation site (https://docs.pangraph.org/), then fast forward `release-docs` to `master`:

```bash
git checkout release-docs
git merge --ff-only origin/master
git push origin release-docs
```

You can combine it with the software release as well (see [docs/dev/developer_guide.md](./docs/dev/developer_guide.md)).

> ️⚠️ Make sure your local branches are up-to-date with the remote repo and that they don't have incompatible changes on them. Otherwise the fast-forward merge will fail. When manipulating `release-*` branches, consider using a git GUI application (e.g. GitKraken) to visually check that you do what you think you do. 
