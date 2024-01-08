ARG DOCKER_BASE_IMAGE
ARG CLANG_VERSION

FROM $DOCKER_BASE_IMAGE as base

SHELL ["bash", "-euxo", "pipefail", "-c"]

ARG DASEL_VERSION="1.22.1"
ARG WATCHEXEC_VERSION="1.17.1"
ARG NODEMON_VERSION="2.0.15"
ARG YARN_VERSION="1.22.18"
ARG CLANG_VERSION=$CLANG_VERSION

RUN set -euxo pipefail >/dev/null \
&& if grep wheezy /etc/apt/sources.list; then export IS_DEBIAN_WHEEZY=1; else export IS_DEBIAN_WHEEZY=0; fi \
&& if [ ${IS_DEBIAN_WHEEZY} == 1 ]; then printf "deb http://archive.debian.org/debian wheezy main non-free contrib\ndeb http://archive.debian.org/debian-security wheezy/updates main non-free contrib\n" > "/etc/apt/sources.list"; fi \
&& if [ ${IS_DEBIAN_WHEEZY} == 1 ]; then echo "Acquire::Check-Valid-Until false;" >> "/etc/apt/apt.conf.d/10-nocheckvalid"; fi \
&& export DEBIAN_FRONTEND=noninteractive \
&& apt-get update -qq --yes \
&& apt-get install -qq --no-install-recommends --yes \
  apt-transport-https \
  bash \
  bash-completion \
  build-essential \
  ca-certificates \
  curl \
  gfortran \
  git \
  gnupg \
  libssl-dev \
  lsb-release \
  make \
  pkg-config \
  python3 \
  python3-pip \
  sudo \
  time \
  xz-utils \
  >/dev/null \
&& \
  if [ "$(lsb_release -cs)" == "wheezy" ]; then \
    echo "deb https://apt.llvm.org/$(lsb_release -cs)/ llvm-toolchain-$(lsb_release -cs) main" >> "/etc/apt/sources.list.d/llvm.list"; \
  else \
    echo "deb https://apt.llvm.org/$(lsb_release -cs)/ llvm-toolchain-$(lsb_release -cs)-${CLANG_VERSION} main" >> "/etc/apt/sources.list.d/llvm.list"; \
  fi \
&& curl -fsSL "https://apt.llvm.org/llvm-snapshot.gpg.key" | sudo apt-key add - \
&& export DEBIAN_FRONTEND=noninteractive \
&& apt-get update -qq --yes \
&& apt-get install -qq --no-install-recommends --yes \
  clang-${CLANG_VERSION} \
  lldb-${CLANG_VERSION} \
  llvm-${CLANG_VERSION} \
  llvm-${CLANG_VERSION}-dev \
  llvm-${CLANG_VERSION}-tools \
  >/dev/null \
&& if [ "$(lsb_release -cs)" != "wheezy" ]; then \
    apt-get install -qq --no-install-recommends --yes \
      clang-tools-${CLANG_VERSION} \
      lld-${CLANG_VERSION} \
      pigz \
      pixz \
    >/dev/null; \
  fi \
&& apt-get clean autoclean >/dev/null \
&& apt-get autoremove --yes >/dev/null \
&& rm -rf /var/lib/apt/lists/*

ARG USER=user
ARG GROUP=user
ARG UID
ARG GID

ENV USER=$USER
ENV GROUP=$GROUP
ENV UID=$UID
ENV GID=$GID
ENV TERM="xterm-256color"
ENV HOME="/home/${USER}"
ENV CARGO_HOME="${HOME}/.cargo"
ENV CARGO_INSTALL_ROOT="${HOME}/.cargo/install"
ENV RUSTUP_HOME="${HOME}/.rustup"
ENV NODE_DIR="/opt/node"
ENV PATH="/usr/lib/llvm-${CLANG_VERSION}/bin:${NODE_DIR}/bin:${HOME}/.local/bin:${HOME}/.cargo/bin:${HOME}/.cargo/install/bin:${PATH}"

# Install Python dependencies
RUN set -euxo pipefail >/dev/null \
&& if [ "$(lsb_release -cs)" == "wheezy" ]; then \
     curl -fsSL https://bootstrap.pypa.io/pip/3.2/get-pip.py | python3 \
   ;fi

RUN set -euxo pipefail >/dev/null \
&& pip3 install --user --upgrade cram

# Install dasel, a tool to query TOML files
RUN set -euxo pipefail >/dev/null \
&& curl -fsSL "https://github.com/TomWright/dasel/releases/download/v${DASEL_VERSION}/dasel_linux_amd64" -o "/usr/bin/dasel" \
&& chmod +x "/usr/bin/dasel" \
&& dasel --version

# Install watchexec - file watcher
RUN set -euxo pipefail >/dev/null \
&& curl -fsSL "https://github.com/watchexec/watchexec/releases/download/cli-v${WATCHEXEC_VERSION}/watchexec-${WATCHEXEC_VERSION}-x86_64-unknown-linux-musl.tar.xz" | tar -C "/usr/bin/" -xJ --strip-components=1 "watchexec-${WATCHEXEC_VERSION}-x86_64-unknown-linux-musl/watchexec" \
&& chmod +x "/usr/bin/watchexec" \
&& watchexec --version

# Install Node.js
COPY .nvmrc /
RUN set -eux >dev/null \
&& if [ "$(lsb_release -cs)" != "wheezy" ]; then \
  mkdir -p "${NODE_DIR}" \
  && cd "${NODE_DIR}" \
  && NODE_VERSION=$(cat /.nvmrc) \
  && curl -fsSL  "https://nodejs.org/dist/v${NODE_VERSION}/node-v${NODE_VERSION}-linux-x64.tar.xz" | tar -xJ --strip-components=1 \
  && npm install -g nodemon@${NODEMON_VERSION} yarn@${YARN_VERSION} >/dev/null \
  && npm config set scripts-prepend-node-path auto \
;fi

# Make a user and group
RUN set -euxo pipefail >/dev/null \
&& \
  if [ -z "$(getent group ${GID})" ]; then \
    groupadd --system --gid ${GID} ${GROUP}; \
  else \
    groupmod -n ${GROUP} $(getent group ${GID} | cut -d: -f1); \
  fi \
&& \
  if [ -z "$(getent passwd ${UID})" ]; then \
    useradd \
      --system \
      --create-home --home-dir ${HOME} \
      --shell /bin/bash \
      --gid ${GROUP} \
      --groups sudo \
      --uid ${UID} \
      ${USER}; \
  fi \
&& sed -i /etc/sudoers -re 's/^%sudo.*/%sudo ALL=(ALL:ALL) NOPASSWD: ALL/g' \
&& sed -i /etc/sudoers -re 's/^root.*/root ALL=(ALL:ALL) NOPASSWD: ALL/g' \
&& sed -i /etc/sudoers -re 's/^#includedir.*/## **Removed the include directive** ##"/g' \
&& echo "foo ALL=(ALL) NOPASSWD: ALL" >> /etc/sudoers \
&& touch ${HOME}/.hushlogin \
&& chown -R ${UID}:${GID} "${HOME}"


USER ${USER}

# Install rustup and toolchain from rust-toolchain.toml
COPY rust-toolchain.toml "${HOME}/rust-toolchain.toml"
RUN set -euxo pipefail >/dev/null \
&& cd "${HOME}" \
&& RUST_TOOLCHAIN=$(dasel select -p toml -s ".toolchain.channel" -f "${HOME}/rust-toolchain.toml") \
&& curl --proto '=https' -sSf https://sh.rustup.rs > rustup-init \
&& chmod +x rustup-init \
&& ./rustup-init -y --no-modify-path --default-toolchain="${RUST_TOOLCHAIN}" \
&& rm rustup-init

# Install toolchain from rust-toolchain.toml and make it default
RUN set -euxo pipefail >/dev/null \
&& cd "${HOME}" \
&& RUST_TOOLCHAIN=$(dasel select -p toml -s ".toolchain.channel" -f "rust-toolchain.toml") \
&& rustup toolchain install "${HOME}" \
&& rustup default "${RUST_TOOLCHAIN}"

# Install remaining toolchain components from rust-toolchain.toml
RUN set -euxo pipefail >/dev/null \
&& cd "${HOME}" \
&& RUST_TOOLCHAIN=$(dasel select -p toml -s ".toolchain.channel" -f "rust-toolchain.toml") \
&& rustup show \
&& rustup default "${RUST_TOOLCHAIN}"


RUN set -euxo pipefail >/dev/null \
&& export SEQKIT_VERSION="2.5.0" \
&& curl -sSL "https://github.com/shenwei356/seqkit/releases/download/v${SEQKIT_VERSION}/seqkit_linux_amd64.tar.gz" | tar -C "${CARGO_HOME}/bin" -xz "seqkit" \
&& chmod +x "${CARGO_HOME}/bin/seqkit"

RUN set -euxo pipefail >/dev/null \
&& export CARGO_BINSTALL_VERSION="1.0.0" \
&& curl -sSL "https://github.com/cargo-bins/cargo-binstall/releases/download/v${CARGO_BINSTALL_VERSION}/cargo-binstall-x86_64-unknown-linux-gnu.tgz" | tar -C "${CARGO_HOME}/bin" -xz "cargo-binstall" \
&& chmod +x "${CARGO_HOME}/bin/cargo-binstall"

RUN set -euxo pipefail >/dev/null \
&& export CARGO_QUICKINSTALL_VERSION="0.2.9" \
&& curl -sSL "https://github.com/alsuren/cargo-quickinstall/releases/download/cargo-quickinstall-${CARGO_QUICKINSTALL_VERSION}-x86_64-unknown-linux-musl/cargo-quickinstall-${CARGO_QUICKINSTALL_VERSION}-x86_64-unknown-linux-musl.tar.gz" | tar -C "${CARGO_HOME}/bin" -xz "cargo-quickinstall" \
&& chmod +x "${CARGO_HOME}/bin/cargo-quickinstall"

RUN set -euxo pipefail >/dev/null \
&& export WASM_BINDGEN_CLI_VERSION="0.2.87" \
&& curl -sSL "https://github.com/rustwasm/wasm-bindgen/releases/download/${WASM_BINDGEN_CLI_VERSION}/wasm-bindgen-${WASM_BINDGEN_CLI_VERSION}-x86_64-unknown-linux-musl.tar.gz" | tar -C "${CARGO_HOME}/bin" --strip-components=1 -xz "wasm-bindgen-${WASM_BINDGEN_CLI_VERSION}-x86_64-unknown-linux-musl/wasm-bindgen" \
&& chmod +x "${CARGO_HOME}/bin/wasm-bindgen"

RUN set -euxo pipefail >/dev/null \
&& export BINARYEN_VERSION="114" \
&& curl -sSL "https://github.com/WebAssembly/binaryen/releases/download/version_${BINARYEN_VERSION}/binaryen-version_${BINARYEN_VERSION}-x86_64-linux.tar.gz" | tar -C "${CARGO_HOME}/bin" --strip-components=2 -xz --wildcards "binaryen-version_${BINARYEN_VERSION}/bin/"'wasm*' \
&& chmod +x ${CARGO_HOME}/bin/wasm*

RUN set -euxo pipefail >/dev/null \
&& export WASM_PACK_VERSION="0.12.1" \
&& curl -sSL "https://github.com/rustwasm/wasm-pack/releases/download/v${WASM_PACK_VERSION}/wasm-pack-v${WASM_PACK_VERSION}-x86_64-unknown-linux-musl.tar.gz" | tar -C "${CARGO_HOME}/bin" --strip-components=1 -xz "wasm-pack-v${WASM_PACK_VERSION}-x86_64-unknown-linux-musl/wasm-pack" \
&& chmod +x "${CARGO_HOME}/bin/wasm-pack"

RUN set -euxo pipefail >/dev/null \
&& export CARGO_WATCH_VERSION="8.4.0" \
&& curl -sSL "https://github.com/watchexec/cargo-watch/releases/download/v${CARGO_WATCH_VERSION}/cargo-watch-v${CARGO_WATCH_VERSION}-x86_64-unknown-linux-gnu.tar.xz" | tar -C "${CARGO_HOME}/bin" --strip-components=1 -xJ "cargo-watch-v${CARGO_WATCH_VERSION}-x86_64-unknown-linux-gnu/cargo-watch" \
&& chmod +x "${CARGO_HOME}/bin/cargo-watch"

USER 0

# Install mmseqs
RUN set -euxo pipefail >/dev/null \
&& export MMSEQS_VERSION="13-45111" \
&& curl -fsSL "https://github.com/soedinglab/MMseqs2/releases/download/${MMSEQS_VERSION}/mmseqs-linux-sse2.tar.gz" | tar -C "/usr/bin/" -xz --strip-components=2 "mmseqs/bin/mmseqs" \
&& mmseqs --help | grep "Version"

# Install mash
RUN set -euxo pipefail >/dev/null \
&& curl -fsSL "https://github.com/marbl/Mash/releases/download/v2.2/mash-Linux64-v2.2.tar" | tar -C "/usr/bin/" -x --strip-components=1 "mash-Linux64-v2.2/mash" \
&& mash --version

# Install fasttree
RUN set -euxo pipefail >/dev/null \
&& curl -fsSLo "/usr/bin/fasttree" "http://www.microbesonline.org/fasttree/FastTree" \
&& chmod +x "/usr/bin/fasttree" \
&& fasttree -help

# Install mafft
RUN set -euxo pipefail >/dev/null \
&& mkdir -p "/opt/mafft" \
&& curl -fsSL "https://mafft.cbrc.jp/alignment/software/mafft-7.490-linux.tgz" \
  | tar xz -C "/usr/bin/" --strip-components=1 "mafft-linux64/" \
&& ln -s "/usr/bin/mafft.bat" "/usr/bin/mafft" \
&& chmod 0777 "/usr/bin/mafft.bat" \
&& chown root:root "/usr/bin/mafft.bat" \
 && chmod -R 0777 "/usr/bin/mafftdir" \
&& chown -R root:root "/usr/bin/mafftdir" \
&& mafft --version

USER ${USER}

# Setup bash
RUN set -euxo pipefail >/dev/null \
&& echo 'alias ll="ls --color=always -alFhp"' >> ~/.bashrc \
&& echo 'alias la="ls -Ah"' >> ~/.bashrc \
&& echo 'alias l="ls -CFh"' >> ~/.bashrc \
&& echo 'function mkcd() { mkdir -p ${1} && cd ${1} ; }' >> ~/.bashrc \
&& rustup completions bash >> ~/.bash_completion \
&& rustup completions bash cargo >>  ~/.bash_completion \
&& echo "source ~/.bash_completion" >> ~/.bashrc

USER ${USER}

WORKDIR ${HOME}/src



# Native compilation for Linux x86_64 with gnu-libc
FROM base as dev

ENV CC_x86_64-unknown-linux-gnu=clang
ENV CXX_x86_64-unknown-linux-gnu=clang++

# Cross-compilation for Linux x86_64 with gnu-libc.
# Same as native, but convenient to have for mass cross-compilation.
FROM dev as cross-x86_64-unknown-linux-gnu

ENV CC_x86_64-unknown-linux-gnu=clang
ENV CXX_x86_64-unknown-linux-gnu=clang++


# Cross-compilation for Linux x86_64 with libmusl
FROM base as cross-x86_64-unknown-linux-musl

USER 0

SHELL ["bash", "-euxo", "pipefail", "-c"]

RUN set -euxo pipefail >/dev/null \
&& curl -fsSL "https://more.musl.cc/11/x86_64-linux-musl/x86_64-linux-musl-cross.tgz" | tar -C "/usr" -xz --strip-components=1

USER ${USER}

RUN set -euxo pipefail >/dev/null \
&& rustup target add x86_64-unknown-linux-musl

ENV CC_x86_64_unknown_linux_musl=x86_64-linux-musl-gcc
ENV CXX_x86_64_unknown_linux_musl=x86_64-linux-musl-g++
ENV CARGO_TARGET_X86_64_UNKNOWN_LINUX_MUSL_LINKER=x86_64-linux-musl-gcc


# Cross-compilation to WebAssembly
FROM base as cross-wasm32-unknown-unknown

SHELL ["bash", "-euxo", "pipefail", "-c"]

RUN set -euxo pipefail >/dev/null \
&& rustup target add wasm32-unknown-unknown


# Cross-compilation for Linux ARM64
FROM base as cross-aarch64-unknown-linux-gnu

USER 0

SHELL ["bash", "-euxo", "pipefail", "-c"]

RUN set -euxo pipefail >/dev/null \
&& export DEBIAN_FRONTEND=noninteractive \
&& apt-get update -qq --yes \
&& apt-get install -qq --no-install-recommends --yes \
  binutils-aarch64-linux-gnu \
  g++-aarch64-linux-gnu \
  gcc-aarch64-linux-gnu \
  gfortran-aarch64-linux-gnu \
  libc6-dev-arm64-cross \
>/dev/null \
&& apt-get clean autoclean >/dev/null \
&& apt-get autoremove --yes >/dev/null \
&& rm -rf /var/lib/apt/lists/*

USER ${USER}

RUN set -euxo pipefail >/dev/null \
&& rustup target add aarch64-unknown-linux-gnu

ENV CC_aarch64_unknown_linux_gnu=aarch64-linux-gnu-gcc
ENV CXX_aarch64_unknown_linux_gnu=aarch64-linux-gnu-g++
ENV CARGO_TARGET_AARCH64_UNKNOWN_LINUX_GNU_LINKER=aarch64-linux-gnu-gcc

# Cross-compilation for Linux ARM64 with libmusl
FROM base as cross-aarch64-unknown-linux-musl

USER 0

SHELL ["bash", "-euxo", "pipefail", "-c"]

RUN set -euxo pipefail >/dev/null \
&& curl -fsSL "https://more.musl.cc/11/x86_64-linux-musl/aarch64-linux-musl-cross.tgz" | tar -C "/usr" -xz --strip-components=1

USER ${USER}

RUN set -euxo pipefail >/dev/null \
&& rustup target add aarch64-unknown-linux-musl

ENV CC_aarch64_unknown_linux_musl=aarch64-linux-musl-gcc
ENV CXX_aarch64_unknown_linux_musl=aarch64-linux-musl-g++
ENV CARGO_TARGET_AARCH64_UNKNOWN_LINUX_MUSL_LINKER=aarch64-linux-musl-gcc


# Cross-compilation for Windows x86_64
FROM base as cross-x86_64-pc-windows-gnu

SHELL ["bash", "-euxo", "pipefail", "-c"]

USER 0

RUN set -euxo pipefail >/dev/null \
&& export DEBIAN_FRONTEND=noninteractive \
&& apt-get update -qq --yes \
&& apt-get install -qq --no-install-recommends --yes \
  gcc-mingw-w64-x86-64 \
  g++-mingw-w64-x86-64 \
  gfortran-mingw-w64-x86-64 \
  binutils-mingw-w64-x86-64  \
>/dev/null \
&& apt-get clean autoclean >/dev/null \
&& apt-get autoremove --yes >/dev/null \
&& rm -rf /var/lib/apt/lists/*

USER ${USER}

RUN set -euxo pipefail >/dev/null \
&& rustup target add x86_64-pc-windows-gnu


# Builds osxcross for Mac cross-compiation
FROM base as osxcross

SHELL ["bash", "-euxo", "pipefail", "-c"]

USER 0

ARG OSXCROSS_URL

# Install cargo-quickinstall
RUN set -euxo pipefail >/dev/null \
&& mkdir -p "/opt/osxcross" \
&& curl -fsSL "${OSXCROSS_URL}" | tar -C "/opt/osxcross" -xJ

USER ${USER}


# Cross-compilation for macOS x86_64
FROM osxcross as cross-x86_64-apple-darwin

SHELL ["bash", "-euxo", "pipefail", "-c"]

USER ${USER}

RUN set -euxo pipefail >/dev/null \
&& rustup target add x86_64-apple-darwin

ENV PATH="/opt/osxcross/bin/:${PATH}"
ENV CC_x86_64-apple-darwin=x86_64-apple-darwin20.2-clang
ENV CXX_x86_64-apple-darwin=x86_64-apple-darwin20.2-clang++
ENV CARGO_TARGET_X86_64_APPLE_DARWIN_LINKER=x86_64-apple-darwin20.2-clang
ENV CARGO_TARGET_X86_64_APPLE_DARWIN_STRIP=x86_64-apple-darwin20.2-strip


# Cross-compilation for macOS ARM64
FROM osxcross as cross-aarch64-apple-darwin

SHELL ["bash", "-euxo", "pipefail", "-c"]

USER ${USER}

RUN set -euxo pipefail >/dev/null \
&& rustup target add aarch64-apple-darwin

ENV PATH="/opt/osxcross/bin/:${PATH}"
ENV CC_aarch64-apple-darwin=aarch64-apple-darwin20.2-clang
ENV CXX_aarch64-apple-darwin=aarch64-apple-darwin20.2-clang++
ENV CARGO_TARGET_AARCH64_APPLE_DARWIN_LINKER=aarch64-apple-darwin20.2-clang
ENV CARGO_TARGET_AARCH64_APPLE_DARWIN_STRIP=aarch64-apple-darwin20.2-strip
