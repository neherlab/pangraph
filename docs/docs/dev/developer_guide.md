---
sidebar_position: 1
---

# Contributing to Pangraph

This guide describe how to setup developer environment, how to build Pangraph, contribute to the codebase and maintain
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

````

## Obtain source code

Pangraph is an open-source project and its source code is available on GitHub under MIT license.

To obtain source code use `git` to clone GitHub repository:

```bash
git clone https://github.com/neherlab/pangraph
```

If you are a team member, use SSH url to be able to push securely and without password prompt
(More details: https://docs.github.com/en/authentication/connecting-to-github-with-ssh)

```bash
git clone git@github.com:neherlab/pangraph.git
```

If you are not a team member, but want to contribute, make a
[fork](https://docs.github.com/en/get-started/quickstart/fork-a-repo), and clone your forked repository instead. You can
then submit changes to pangraph as
a [Pull Request](https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/creating-a-pull-request)
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

We run tests using `cargo-nextest` (https://nexte.st/). You can install it from [GitHub Releases](https://github.com/nextest-rs/nextest/releases) or build and install from source with

```bash
cargo install cargo-nextest --locked
```

#### All tests

Run all tests with:

```bash
cargo nextest run
```

Add `--no-fail-fast` flag to keep going even if there are failures.

A subset of tests can be ran by providing a regex matching full test name. For example to run tests `test_foo` and `test_foo_bar`

```bash
cargo nextest run foo
```

You may experiment with command line arguments and the `prettytest` script for the most useful and pleasant output:

```bash
cargo -q nextest run --success-output=immediate --workspace --cargo-quiet --no-fail-fast --hide-progress-bar --color=always | "./dev/prettytest"
```

Arguments of `cargo test` (and by extension `cargo nextest`) are somewhat confusing, so you could setup some personal shell scripts or aliases to simplify routine work on Rust projects.

See also: 
 
 - [cargo-test docs](https://doc.rust-lang.org/cargo/commands/cargo-test.html)
 - [cargo-nextest docs](https://nexte.st/)
 - Read more about different test types in [The Rust Book](https://doc.rust-lang.org/book/Rch11-03-test-organization.html) as well as in [Rust by Example](https://doc.rust-lang.org/rust-by-example/testing.html).

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

Clippy is configured in `clippy.toml` and `.cargo/config.toml`.

### Formatting (code style)

Code formatting is done using [rustfmt](https://rust-lang.github.io/rustfmt):

```bash
cargo fmt --all
```

Rustfmt is configured in `rustfmt.toml`.

## Maintenance

### Upgrading Rust

Rust version is defined in `rust-toolchain.toml`. When using `cargo`, the version defined in this file gets installed
automatically.

### Upgrading Rust

Dependencies for subprojects are defined in  `packages/**/Cargo.toml` and in `Cargo.lock`. They are periodically
upgraded by a dedicated maintainer, manually using `cargo-upgrade`
from [cargo-edit](https://github.com/killercup/cargo-edit) package.

```bash
cargo upgrade --workspace
```

### Documentation

End-user and developer documentation is in `docs/` subdirectory. For details, see `README.md` file there.


### Versioning

TODO

### Releases

TODO

### Continuous integration (CI), packaging and distribution

TODO
