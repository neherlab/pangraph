# Stage: builder image
# This istage builds use a lot of dependencies and produce the binaries. The results will be copied
# to another image and the builder image will be discarded.
FROM debian:11 as builder

SHELL ["bash", "-c"]


RUN set -euxo pipefail \
&& export DEBIAN_FRONTEND=noninteractive \
&& apt-get update -qq --yes \
&& apt-get install -qq --no-install-recommends --yes \
  build-essential \
  ca-certificates \
  curl \
  mafft \
  make \
  mash \
>/dev/null \
&& apt-get autoremove --yes >/dev/null \
&& apt-get clean autoclean >/dev/null \
&& rm -rf /var/lib/apt/lists/*

ENV PATH="/build_dir/bin:/build_dir/vendor/julia/bin:$PATH"

COPY . /build_dir/

RUN set -euxo pipefail \
&& cd /build_dir \
&& make


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

# Allows non-root users to read dependencies
RUN set -euxo pipefail \
&& chmod -R +r /root/ \
&& chmod +x /root/

CMD ["/usr/bin/pangraph"]
