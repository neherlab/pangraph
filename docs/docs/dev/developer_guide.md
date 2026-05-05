---
sidebar_position: 1
---

# Contributing to Pangraph

This guide describes how to setup developer environment, how to build Pangraph, contribute to the codebase and maintain
the project.

## Setup developer environment

This guide assumes Ubuntu 24.04 operating system, but will likely work similarly to any other Linux and Unix-like
machine.

Pangraph is written in Rust. The usual `rustup` & `cargo` workflow can be used:

```bash
# Install required dependencies
# These particular commands are specific for Ubuntu Linux and will work on some other Debian-based Linux distros.
# Refer to documentation of your operating system to find how to install these dependencies.
sudo apt-get update
sudo apt-get install \
  bash \
  build-essential \
  clang \
  curl \
  gcc \
  git \
  libbz2-dev \
  libclang-dev \
  liblzma-dev \
  libssl-dev \
  libzstd-dev \
  make \
  pkg-config \
  zlib1g-dev \

# Install Rustup, the Rust version manager (https://www.rust-lang.org/tools/install)
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | bash -s -- -y

# Add Rust tools to the $PATH. You might add this line to your .bashrc or .zshrc so that the $PATH is adjusted automatically when a new terminal session is opened.
export PATH="$PATH:$HOME/.cargo/bin"

```

## Obtain source code

Pangraph is an open-source project and its source code is available on GitHub under MIT license.

To obtain source code use `git` to clone GitHub repository:

```bash
git clone https://github.com/neherlab/pangraph
```

If you are a team member, use SSH url to be able to push securely and without password prompt (More details:
https://docs.github.com/en/authentication/connecting-to-github-with-ssh)

```bash
git clone git@github.com:neherlab/pangraph.git
```

If you are not a team member, but want to contribute, make a
[fork](https://docs.github.com/en/get-started/quickstart/fork-a-repo), and clone your forked repository instead. You can
then submit changes to pangraph as a
[Pull Request](https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/creating-a-pull-request)
. Pangraph maintainers will then review your changes and will consider merging them into the main project.

## Build and run

### Build and run executables

```bash
# Go to the cloned directory
cd pangraph

# (optional) checkout a branch (different from default)
git checkout <branch_name>

# Build and run in debug mode (convenient for development, fast to build, slow to run, has more information in stack traces and when running under a debugger)
cargo run --bin=pangraph

# Run with additional arguments passed to the executable
cargo run --bin=pangraph -- --help

# Instead of `--bin=pangraph` you can also run any other executable from `packages/pangraph/src/bin/`. Just substitute its filename.
# This is a `cargo` convention: everything in `src/bin/` that has a `main()` function in it becomes an executable. This way you can add more executables.
cargo run --bin=my_executable

# Run Pangraph in release mode (slow to build, fast to run, very little information in stack traces and during debugging)
cargo run --release --bin=pangraph

# Alternatively, build and run separately. The compiled binaries will be in `target/` directory by default.
cargo run --release --bin=pangraph
./target/release/pangraph

```

Note, on first build of a particular project, `cargo` will search for one of the possible toolchain config files and
will automatically install Rust version required by the project. This may cause first build to take longer than usual.

### Testing

#### Install requirements

We run tests using `cargo-nextest` (https://nexte.st/). You can install it from
[GitHub Releases](https://github.com/nextest-rs/nextest/releases) or build and install from source with

```bash
cargo install cargo-nextest --locked
```

#### All tests

Run all tests with:

```bash
cargo nextest run
```

Add `--no-fail-fast` flag to keep going even if there are failures.

A subset of tests can be ran by providing a regex matching full test name. For example to run tests `test_foo` and
`test_foo_bar`

```bash
cargo nextest run foo
```

You may experiment with command line arguments and the `prettytest` script for the most useful and pleasant output:

```bash
cargo -q nextest run --success-output=immediate --workspace --cargo-quiet --no-fail-fast --hide-progress-bar --color=always | "./dev/prettytest"
```

Arguments of `cargo test` (and by extension `cargo nextest`) are somewhat confusing, so you could setup some personal
shell scripts or aliases to simplify routine work on Rust projects.

See also:

- [cargo-test docs](https://doc.rust-lang.org/cargo/commands/cargo-test.html)
- [cargo-nextest docs](https://nexte.st/)
- Read more about different test types in [The Rust Book](https://doc.rust-lang.org/book/ch11-03-test-organization.html)
  as well as in [Rust by Example](https://doc.rust-lang.org/rust-by-example/testing.html).

#### Unit tests

If you want to run only unit tests:

```bash
cargo nextest run --lib
cargo nextest run --lib foo
```

#### Integration tests

If you want to run only integration tests:

```bash
cargo nextest run --test='*'
cargo nextest run --test='*' foo
cargo nextest run --test='foo'
```

### Linting (static analysis)

Rust code is linted by running [Clippy](https://github.com/rust-lang/rust-clippy):

```bash
cargo clippy --all-targets --all
```

Clippy is configured in [`.cargo/config.toml`](https://github.com/neherlab/pangraph/blob/master/.cargo/config.toml).

### Formatting (code style)

Code formatting is done using [rustfmt](https://rust-lang.github.io/rustfmt):

```bash
cargo fmt --all
```

Rustfmt is configured in [`rustfmt.toml`](https://github.com/neherlab/pangraph/blob/master/rustfmt.toml).

### Development scripts (optional)

The project includes development scripts at
[`./dev/docker/run`](https://github.com/neherlab/pangraph/blob/master/dev/docker/run) and
[`./dev/dev`](https://github.com/neherlab/pangraph/blob/master/dev/dev) which provide shortcuts for common development
tasks. All commands run inside a Docker container:

```bash
# Build in debug mode
./dev/docker/run ./dev/dev b

# Build in release mode
./dev/docker/run ./dev/dev br

# Run in debug mode
./dev/docker/run ./dev/dev r pangraph -- build --help

# Run in release mode
./dev/docker/run ./dev/dev rr pangraph -- build --help

# Run tests
./dev/docker/run ./dev/dev t

# Run unit tests only
./dev/docker/run ./dev/dev tu

# Run integration tests only
./dev/docker/run ./dev/dev ti

# Run linter
./dev/docker/run ./dev/dev l

# Run linter with auto-fixes
./dev/docker/run ./dev/dev lf

# Format code
./dev/docker/run ./dev/dev f

# Run arbitrary command in the container
./dev/docker/run your command here
```

The same docker commands and the same container is used in CI. This setup ensures a consistent, isolated, reproducible
environment on local machines and on remotes.

Additional cargo aliases are defined in
[`.cargo/config.toml`](https://github.com/neherlab/pangraph/blob/master/.cargo/config.toml) and can be used directly
with `cargo` (e.g., `cargo l` for lint, `cargo t` for tests). These are optional shortcuts - canonical `cargo` commands
work as usual.

## Maintenance

### Upgrading Rust

Rust version is defined in
[`rust-toolchain.toml`](https://github.com/neherlab/pangraph/blob/master/rust-toolchain.toml). When using `cargo`, the
version defined in this file gets installed automatically.

### Upgrading Rust dependencies

Dependencies for subprojects are defined in `packages/**/Cargo.toml` and in
[`Cargo.lock`](https://github.com/neherlab/pangraph/blob/master/Cargo.lock). They are periodically upgraded by a
dedicated maintainer, manually using `cargo-upgrade` from [cargo-edit](https://github.com/killercup/cargo-edit) package.

```bash
cargo upgrade --workspace
```

or, for a more radical upgrade:

```bash
cargo upgrade --pinned --incompatible --verbose --recursive
```

Then remove `Cargo.lock` and rebuild the project to apply the upgrades across dependency tree.

### Documentation

End-user and developer documentation source is in [`docs/`](../../), built using [Docusaurus](https://docusaurus.io/).

#### Requirements

- Node.js via [nvm](https://github.com/nvm-sh/nvm). The supported version is declared in
  [`.nvmrc`](https://github.com/neherlab/pangraph/blob/master/.nvmrc).
- [bun](https://bun.sh)

#### Local development

```bash
cd docs

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

The development server serves the documentation website on local port `4000`. Edits to content and code are reflected on
the fly.

Documentation files are written in markdown (`.md` or `.mdx`) and located under [`docs/docs/`](../).

#### Production build

To build and inspect the production build:

```bash
cd docs
bun prod:build
bun prod:serve
```

The `build` command produces HTML, CSS and JS files ready for deployment. The `serve` command serves them on port
`5000`.

#### Generate command-line reference

The [`docs/docs/reference.md`](../reference.md) file is generated using script
[`generate-reference-docs`](https://github.com/neherlab/pangraph/blob/master/docs/generate-reference-docs):

```bash
cargo build --bin=pangraph
cd docs
./generate-reference-docs "../target/debug/pangraph" "docs/reference.md"
```

Do not edit the generated file manually. All manual changes will be overwritten by automation.

### Versioning and releases

There are multiple release targets. Each has a dedicated script in `dev/` and a corresponding release branch that
triggers CI on push.

| Target     | Script                     | Branch               | CI workflow                                                                                         | Destination                                                                                                               |
| ---------- | -------------------------- | -------------------- | --------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------- |
| CLI        | `./dev/release`            | `release-cli`        | [cli.yml](https://github.com/neherlab/pangraph/blob/master/.github/workflows/cli.yml)               | [GitHub Releases](https://github.com/neherlab/pangraph/releases), [DockerHub](https://hub.docker.com/r/neherlab/pangraph) |
| PyPangraph | `./dev/release-pypangraph` | `release-pypangraph` | [pypangraph.yml](https://github.com/neherlab/pangraph/blob/master/.github/workflows/pypangraph.yml) | [PyPI](https://pypi.org/project/pypangraph/)                                                                              |
| Docs       | `./dev/release-docs`       | `release-docs`       | [docs.yml](https://github.com/neherlab/pangraph/blob/master/.github/workflows/docs.yml)             | [docs.pangraph.org](https://docs.pangraph.org)                                                                            |

#### Releasing Pangraph CLI

1. Check out `master` branch. Make sure that all the changes that you want to release are on `master` branch. Make sure
   that the last [GitHub Action](https://github.com/neherlab/pangraph/actions) on master branch succeeded (use branch
   filter dropdown). If not, make sure to fix the failures before trying to release.

2. Prepare changelog document for the release: open
   [`CHANGELOG.md`](https://github.com/neherlab/pangraph/blob/master/CHANGELOG.md) in the root directory of the project,
   add `## Unreleased` section at the top of the file, spelled exactly like this - important for automation. Under this
   section, describe all changes in the coming release. This is a user-facing document, so use simple words, avoid
   internal and dev jargon. If `## Unreleased` section already exists, then extend it - do not add multiple of these
   sections.

3. Perform pre-release checks, bump versions, commit using the helper
   [`./dev/release`](https://github.com/neherlab/pangraph/blob/master/dev/release) script.

   Read comments in the script on how to install currently required dependencies.

   This step is local to your working copy of the project, it does not push or otherwise publishes anything yet.

   ```bash
   ./dev/release ${bump_type}
   ```

   most likely you want "bump type" to be one of :
   - `patch` - for bug fix releases, e.g. `1.1.0` -> `1.1.1`
   - `minor` - for new feature releases, e.g. `1.1.0` -> `1.2.0`
   - `major` - for releases with breaking changes, e.g. `1.1.0` -> `2.0.0`

   Having hard times to decide? Read [Semantic Versioning](https://semver.org/).

4. Follow instructions printed by the script. Resolve errors, if any. If finished successfully, follow instructions on
   how to fast-forward and push the changes to the `release-cli` branch.

5. Optionally, you can combine the releases of CLI, docs and PyPangraph. Just add more commits to `master` branch and
   then fast-forward and push them to corresponding branches all together.

6. The push to `release-cli` triggers a GitHub Action. This is non-reversible step. Do the push. Watch for any
   malfunctions in [GitHub Actions](https://github.com/neherlab/pangraph/actions).

7. If GitHub Action succeeds, double check that release is up on
   [GitHub Releases](https://github.com/neherlab/pangraph/releases) and
   [Docker Hub](https://hub.docker.com/r/neherlab/pangraph).

#### Releasing PyPangraph

1. Check out `master` branch. Make sure that all the changes you want to release are on `master` branch.

2. Prepare changelog for the release: open
   [`packages/pypangraph/CHANGELOG.md`](https://github.com/neherlab/pangraph/blob/master/packages/pypangraph/CHANGELOG.md)
   and add a `## ${version_number}` section describing released changes. The changelog does not need to be committed -
   the release script will include it in the release commit.

3. Run the release script:

   ```bash
   ./dev/release-pypangraph ${version_number}
   ```

   The script validates preconditions, bumps the version in
   [`pyproject.toml`](https://github.com/neherlab/pangraph/blob/master/packages/pypangraph/pyproject.toml), commits both
   the version bump and changelog together, and fast-forwards `release-pypangraph` to HEAD.

   Dependencies: `dasel` (see comments in
   [`./dev/release`](https://github.com/neherlab/pangraph/blob/master/dev/release) for install instructions).

4. Follow instructions printed by the script. Push master and the release branch.

5. The push to `release-pypangraph` triggers a GitHub Action. This is a non-reversible step. Watch for any malfunctions
   in [GitHub Actions](https://github.com/neherlab/pangraph/actions).

6. If GitHub Action succeeds, verify that the release is up on [PyPI](https://pypi.org/project/pypangraph).

#### Releasing user documentation

Documentation is deployed to [docs.pangraph.org](https://docs.pangraph.org) by pushing the `release-docs` branch, which
triggers the [docs.yml](https://github.com/neherlab/pangraph/blob/master/.github/workflows/docs.yml) workflow.

1. Make sure documentation changes are on `master` branch and CI passes.

2. Run the release script:

   ```bash
   ./dev/release-docs
   ```

   The script fetches origin, verifies there are documentation changes since the last release, and fast-forwards
   `release-docs` to `origin/master`.

3. Follow instructions printed by the script. Push the release branch.

4. The push triggers CI which builds and deploys to AWS S3/CloudFront and creates a timestamped git tag
   (`docs-YYYY-MM-DD_HH-MM-SSZ`). Watch for malfunctions in
   [GitHub Actions](https://github.com/neherlab/pangraph/actions).

5. Verify the changes are live at [docs.pangraph.org](https://docs.pangraph.org).
