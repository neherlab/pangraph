# Stage: builder image
# This istage builds use a lot of dependencies and produce the binaries. The results will be copied
# to another image and the builder image will be discarded.
FROM ubuntu:20.04 as builder

SHELL ["bash", "-c"]


RUN set -euxo pipefail \
&& export DEBIAN_FRONTEND=noninteractive \
&& apt-get update -qq --yes \
&& apt-get install -qq --no-install-recommends --yes \
  build-essential \
  ca-certificates \
  curl \
  make \
  mafft \
>/dev/null \
&& apt-get autoremove --yes >/dev/null \
&& apt-get clean autoclean >/dev/null \
&& rm -rf /var/lib/apt/lists/*

# TODO: We need to set the PATH to Julia bin dir. However the version is hardwired into the path.
# We need to install Julia to a version-neutral dir.
ENV PATH="/build_dir/bin:/build_dir/vendor/julia/bin:$PATH"

COPY bin /build_dir/bin

RUN set -euxo pipefail \
&& mkdir -p /build_dir/vendor

COPY . /build_dir/

RUN set -euxo pipefail \
&& cd /build_dir \
&& jc=$(which julia) make


# Stage: production image
# We start over, from clean debian image, and copy the binaries from the builder stage.
FROM debian:11 as prod

# Copy pangraph from the builder stage
COPY --from=builder /build_dir/pangraph/ /usr/

# Copy julia dependencies from the builder stage
COPY --from=builder /root/.julia/artifacts /root/.julia/artifacts
COPY --from=builder /root/.julia/conda/3/bin /root/.julia/conda/3/bin
COPY --from=builder /root/.julia/conda/3/lib /root/.julia/conda/3/lib

SHELL ["bash", "-c"]

RUN set -euxo pipefail \
&& export DEBIAN_FRONTEND=noninteractive \
&& apt-get update -qq --yes \
&& apt-get install -qq --no-install-recommends --yes \
  mafft \
  mash \
>/dev/null \
&& apt-get autoremove --yes >/dev/null \
&& apt-get clean autoclean >/dev/null \
&& rm -rf /var/lib/apt/lists/*


CMD ["/usr/bin/pangraph"]
