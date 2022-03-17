#FROM julia:1.7.2-bullseye as dev
FROM ubuntu:20.04 as builder

SHELL ["bash", "-c"]


RUN set -euxo pipefail \
&& export DEBIAN_FRONTEND=noninteractive \
&& apt-get update -qq --yes \
&& apt-get install -qq --no-install-recommends --yes \
  bash \
  build-essential \
  ca-certificates \
  coreutils \
  curl \
  git \
  make \
  python3 \
  python3-dev \
  python3-pip \
  python3-setuptools \
>/dev/null \
&& apt-get autoremove --yes >/dev/null \
&& apt-get clean autoclean >/dev/null \
&& rm -rf /var/lib/apt/lists/*

# TODO: We need to set the PATH to Julia bin dir. However the version is hardwired into the path.
# We need to install Julia to a version-neutral dir.
ENV PATH="/build_dir/bin:/build_dir/vendor/julia-1.7.2/bin:$PATH"

ENV JULIA_DEPOT_PATH="/build_dir/.cache/julia"

COPY . /build_dir/

RUN set -euxo pipefail \
&& cd /build_dir \
&& ./bin/setup-pangraph

# HACK: Cannot call `make` at this point, because it will not call `Pkg.build()` and `Pkg.instantiate()`.
# So downloading Julia here temporarily to be able to call them.
RUN set -euxo pipefail \
&& cd /build_dir/vendor \
&& curl -L https://julialang-s3.julialang.org/bin/linux/x64/1.7/julia-1.7.2-linux-x86_64.tar.gz -o julia-1.7.2-linux-x86_64.tar.gz \
&& tar xzf julia-1.7.2-linux-x86_64.tar.gz

# TODO: This should be in the Makefile
RUN set -euxo pipefail \
&& cd /build_dir \
&& julia --project=. -e 'import Pkg; Pkg.instantiate()'

# TODO: This should be in the Makefile
RUN set -euxo pipefail \
&& cd /build_dir \
&& julia --project=. -e 'import Pkg; Pkg.build()'

RUN set -euxo pipefail \
&& cd /build_dir \
&& make


FROM debian:11 as prod

COPY --from=builder /build_dir/pangraph/ /usr/

CMD ["/usr/bin/pangraph"]
