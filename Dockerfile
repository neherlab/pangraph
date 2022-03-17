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
>/dev/null \
&& apt-get autoremove --yes >/dev/null \
&& apt-get clean autoclean >/dev/null \
&& rm -rf /var/lib/apt/lists/*

# TODO: We need to set the PATH to Julia bin dir. However the version is hardwired into the path.
# We need to install Julia to a version-neutral dir.
ENV PATH="/build_dir/bin:/build_dir/vendor/julia-1.7.2/bin:$PATH"

COPY bin /build_dir/bin

RUN set -euxo pipefail \
&& cd /build_dir \
&& ./bin/setup-pangraph

# HACK: Cannot call `make` at this point, because it will not call `Pkg.build()` and `Pkg.instantiate()`.
# So downloading Julia here temporarily to be able to call them.
RUN set -euxo pipefail \
&& mkdir -p /build_dir/vendor \
&& cd /build_dir/vendor \
&& curl -L https://julialang-s3.julialang.org/bin/linux/x64/1.7/julia-1.7.2-linux-x86_64.tar.gz -o julia-1.7.2-linux-x86_64.tar.gz \
&& tar xzf julia-1.7.2-linux-x86_64.tar.gz

COPY . /build_dir/

# TODO: This should be in the Makefile
RUN set -euxo pipefail \
&& cd /build_dir \
&& julia --project=. -e 'import Pkg; Pkg.instantiate()'

# TODO: This should be in the Makefile
RUN set -euxo pipefail \
&& cd /build_dir \
&& julia --project=. -e 'import Pkg; Pkg.add(name="Conda")' \
&& julia --project=. -e 'import Conda; Conda.add("ete3", channel="etetoolkit")'

# TODO: This should be in the Makefile
RUN set -euxo pipefail \
&& cd /build_dir \
&& julia --project=. -e 'import Pkg; Pkg.build()'

RUN set -euxo pipefail \
&& cd /build_dir \
&& jc=$(which julia) make


# Stage: production image
# We start over, from clean debian image, and copy the binaries from the builder stage.
FROM debian:11 as prod

# Copy dependencies
COPY --from=builder /build_dir/pangraph/ /usr/
COPY --from=builder /root/.julia/conda/3 /root/.julia/conda/3

CMD ["/usr/bin/pangraph"]
