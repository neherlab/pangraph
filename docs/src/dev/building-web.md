# Building Pangraph web app

## Install dependencies

Install Node.js version 18+ (the latest LTS release is recommended), by either downloading it from the official website: https://nodejs.org/en/download/, or by using [nvm](https://github.com/nvm-sh/nvm). We don't recommend using Node.js from the package manager of your operating system, and neither from conda nor other sources. Make sure the `node`, `npm` and `yarn` executables are in your `$PATH`:

```bash
node --version
npm --version
yarn --version
```

## Build and run in development mode

Development mode is suitable only for local day-to-day development. It is relatively fast to build and rebuild (incrementally), slow to run, contains lots of debug info. It looks and feels mostly like the final app the users will see, but the way it is implemented is very different. So dev mode is only suitable for day-to-day dev tasks. All testing must be performed on production app.

In order to build and run the app in dev mode. Run this sequence of commands:

```bash
cd web

# Prepare environment variables with default values
cp .env.example .env

# Install dependency packages
yarn install

# Run development server
yarn dev
```

Open http://localhost:3000 in the browser.

The build is lazy, so page code is only compiled when you make a request that page in the browser. Code modifications should trigger incremental rebuild and fast refresh in the browser.


## Build and run in production mode

Production version of the app closely corresponds to what will be shipped to end users. The build is always a full build and the result of it are the HTML, CSS and JS files which can be served using any static web server.

In order to build and run the app in production mode (fast to build, slow to run, lots of debug info) Run this sequence of commands:

```bash
cd web

# Prepare environment variables with default values
cp .env.example .env

# Install dependency packages
yarn install

# Build production app files
yarn prod:build

# Serve production app files locally
yarn prod:serve

# There is also a shortcut that runs prod:build && prod:serve:
# yarn prod
```

Open a http://localhost:8080 in the browser.

The production version has no fast refresh and it always performs the full rebuild.


## Deployment

TODO


## Feedback

If something does not work, or if you have ideas on improving the build setup, feel free to open an issue or a pull request.
