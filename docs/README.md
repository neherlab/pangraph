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

## Deployment

TODO
