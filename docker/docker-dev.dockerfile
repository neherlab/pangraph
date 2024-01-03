ARG DOCKER_BASE_IMAGE

FROM $DOCKER_BASE_IMAGE as dev

SHELL ["bash", "-euxo", "pipefail", "-c"]

ARG DOCKER_BASE_IMAGE

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
ENV PATH="${NODE_DIR}/bin:${HOME}/.local/bin:${PATH}"

# Install jq, a tool to query JSON files
RUN set -euxo pipefail >/dev/null \
&& curl -fsSL -o "/usr/bin/jq" "https://github.com/stedolan/jq/releases/download/jq-1.6/jq-linux64" \
&& chmod +x "/usr/bin/jq" \
&& jq --version

# Install dasel, a tool to query TOML files
RUN set -euxo pipefail >/dev/null \
&& curl -fsSL "https://github.com/TomWright/dasel/releases/download/v1.22.1/dasel_linux_amd64" -o "/usr/bin/dasel" \
&& chmod +x "/usr/bin/dasel" \
&& dasel --version

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
