ARG DOCKER_BASE_IMAGE

FROM $DOCKER_BASE_IMAGE as dev

SHELL ["bash", "-euxo", "pipefail", "-c"]

ARG DOCKER_BASE_IMAGE
ARG DASEL_VERSION="1.22.1"

# Install required packages if running Debian or Ubuntu
RUN set -euxo pipefail >/dev/null \
&& export DEBIAN_FRONTEND=noninteractive \
&& apt-get update -qq --yes \
&& apt-get install -qq --no-install-recommends --yes \
  apt-transport-https \
  bash \
  bash-completion \
  build-essential \
  ca-certificates \
  curl \
  git \
  gnupg \
  libssl-dev \
  lsb-release \
  parallel \
  pigz \
  pixz \
  pkg-config \
  python3 \
  python3-pip \
  rename \
  sudo \
  time \
  xz-utils \
>/dev/null \
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
ENV NODE_DIR="/opt/node"
ENV PATH="${NODE_DIR}/bin:${HOME}/.local/bin:${HOME}/.cargo/bin:${HOME}/.cargo/install/bin:${PATH}"

# Install dasel, a tool to query TOML files
RUN set -euxo pipefail >/dev/null \
&& curl -fsSL -o "/usr/bin/jq" "https://github.com/stedolan/jq/releases/download/jq-1.6/jq-linux64" \
&& chmod +x "/usr/bin/jq" \
&& jq --version

# Install dasel, a tool to query TOML files
RUN set -euxo pipefail >/dev/null \
&& curl -fsSL "https://github.com/TomWright/dasel/releases/download/v${DASEL_VERSION}/dasel_linux_amd64" -o "/usr/bin/dasel" \
&& chmod +x "/usr/bin/dasel" \
&& dasel --version


COPY docker/files /files

COPY Manifest.toml "${HOME}/Manifest.toml"

RUN set -euxo pipefail >/dev/null \
&& cd "${HOME}" \
&& JULIA_VERSION=$(dasel select -p toml -s ".julia_version" -f "${HOME}/Manifest.toml") \
&& JULIA_VERSION_SHORT=${JULIA_VERSION%.*} \
&& curl -fsSL "https://julialang-s3.julialang.org/bin/linux/x64/${JULIA_VERSION_SHORT}/julia-${JULIA_VERSION}-linux-x86_64.tar.gz" | tar -C "/usr" --strip-components=1 -xz \
&& julia --version

RUN set -euxo pipefail >/dev/null \
&& julia --project="" -e 'import Pkg; Pkg.update(); Pkg.add("IJulia"); Pkg.precompile(); using IJulia; installkernel("Julia", "--project=@.")'

RUN set -euxo pipefail >/dev/null \
&& pip3 install \
  jupyter

# Make a user and group
RUN set -euxo pipefail >/dev/null \
&& \
  if [ -z "$(getent group ${GID})" ]; then \
    groupadd --system --gid ${GID} ${GROUP}; \
  else \
    groupmod -n ${GROUP} $(getent group ${GID} | cut -d: -f1); \
  fi \
&& export SUDO_GROUP="sudo" \
&& \
  if [[ "$DOCKER_BASE_IMAGE" == centos* ]] || [[ "$DOCKER_BASE_IMAGE" == *manylinux2014* ]]; then \
    export SUDO_GROUP="wheel"; \
  fi \
&& \
  if [ -z "$(getent passwd ${UID})" ]; then \
    useradd \
      --system \
      --create-home --home-dir ${HOME} \
      --shell /bin/bash \
      --gid ${GROUP} \
      --groups ${SUDO_GROUP} \
      --uid ${UID} \
      ${USER}; \
  fi \
&& sed -i /etc/sudoers -re 's/^%sudo.*/%sudo ALL=(ALL:ALL) NOPASSWD:ALL/g' \
&& sed -i /etc/sudoers -re 's/^root.*/root ALL=(ALL:ALL) NOPASSWD:ALL/g' \
&& sed -i /etc/sudoers -re 's/^#includedir.*/## **Removed the include directive** ##"/g' \
&& echo "%sudo ALL=(ALL) NOPASSWD:ALL" >> /etc/sudoers \
&& echo "${USER} ALL=(ALL) NOPASSWD:ALL" >> /etc/sudoers \
&& touch ${HOME}/.hushlogin \
&& chown -R ${UID}:${GID} "${HOME}"


USER ${USER}

# Setup bash
RUN set -euxo pipefail >/dev/null \
&& echo 'alias ll="ls --color=always -alFhp"' >> ~/.bashrc \
&& echo 'alias la="ls -Ah"' >> ~/.bashrc \
&& echo 'alias l="ls -CFh"' >> ~/.bashrc \
&& echo 'function mkcd() { mkdir -p ${1} && cd ${1} ; }' >> ~/.bashrc \
&& echo "source ~/.bash_completion" >> ~/.bashrc

USER ${USER}
