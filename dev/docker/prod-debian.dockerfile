# syntax=docker/dockerfile:1
# check=experimental=all
FROM debian:12

SHELL ["bash", "-euxo", "pipefail", "-c"]


RUN set -euxo pipefail \
&& ln -s /usr/bin/pangraph /pangraph \
&& export DEBIAN_FRONTEND=noninteractive \
&& apt-get update -qq --yes \
&& apt-get install -qq --no-install-recommends --yes \
  bash \
  ca-certificates \
  curl \
  tar \
>/dev/null \
&& apt-get autoremove --yes >/dev/null \
&& apt-get clean autoclean >/dev/null \
&& rm -rf /var/lib/apt/lists/*


RUN set -euxo pipefail \
&& export MMSEQS_VERSION="17-b804f" \
&& curl -fsSL "https://github.com/soedinglab/MMseqs2/releases/download/${MMSEQS_VERSION}/mmseqs-linux-sse41.tar.gz" \
  | tar --strip-components=2 -C "/usr/bin" -xz "mmseqs/bin/mmseqs" \
&& ls "/usr/bin/mmseqs" \
&& /usr/bin/mmseqs --help \
&& mmseqs --help


COPY .out/pangraph-x86_64-unknown-linux-gnu /usr/bin/pangraph
RUN set -euxo pipefail \
&& ls "/usr/bin/pangraph" \
&& /usr/bin/pangraph --help \
&& pangraph --help \
