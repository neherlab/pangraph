.PHONY:	all install pangraph environment release documentation clean
.SUFFIXES:
.SECONDARY:

version := 1.7.2

ifeq ($(jc),)
jc := ./vendor/julia/bin/julia
endif

jflags := --project=.
srcs   := $(wildcard src/*.jl src/*/*.jl)
# julia  := julia $(jflags)

datadir   := data/synthetic
testdatum := $(datadir)/test.fa

all: pangraph

install: pangraph/bin/pangraph
	ln -s $$(pwd)/$< bin/pangraph

environment: $(jc)
	$(jc) $(jflags) -e 'using Pkg; pkg"add TreeTools#functions_for_pangraph"; Pkg.instantiate();'
	$(jc) $(jflags) -e 'import Pkg; Pkg.build();'

pangraph: pangraph/bin/pangraph

$(datadir):
	mkdir -p $@

$(testdatum): | environment $(jc) $(datadir)
	$(jc) $(jflags) -e 'using PanGraph; PanGraph.Simulation.test()'

# TODO: look for ARM vs x86
# TODO: julia gets installed into a directory containing version number. This makes it impossible to refer to the
#  installation outside of this file.
$(jc):
ifeq ($(shell uname -s),Linux)
	mkdir -p vendor && \
	cd vendor && \
	curl -L https://julialang-s3.julialang.org/bin/linux/x64/$(basename $(version))/julia-$(version)-linux-x86_64.tar.gz -o julia-$(version)-linux-x86_64.tar.gz && \
	tar xzf julia-$(version)-linux-x86_64.tar.gz && \
	mv julia-$(version) julia
else
ifeq ($(shell uname -s),Darwin)
	mkdir -p vendor && \
	cd vendor && \
	curl -L https://julialang-s3.julialang.org/bin/mac/x64/$(basename $(version))/julia-$(version)-mac64.tar.gz -o julia-$(version)-mac64.tar.gz && \
	tar xzf julia-$(version)-mac64.tar.gz && \
	mv julia-$(version) julia
else
	$(error unsupported host system)
endif
endif

pangraph/bin/pangraph: compile.jl trace.jl $(srcs) $(testdatum) $(jc)
	$(jc) $(jflags) $<

documentation: environment
	$(jc) $(jflags) -e 'import Pkg; Pkg.add(name="Documenter");'
	$(jc) $(jflags) -e 'import Pkg; Pkg.build();'
	$(jc) $(jflags) docs/make.jl

release:
	tar czf pangraph.tar.gz pangraph

clean:
	rm -rf pangraph pangraph.tar.gz


export CONTAINER_NAME=neherlab/pangraph

SHELL=bash
.ONESHELL:
docker:
	set -euxo pipefail

	export DOCKER_TAGS="--tag $${CONTAINER_NAME}:latest"

	if [ ! -z "$${GIT_TAG:-}" ]; then
		export DOCKER_TAGS="$${DOCKER_TAGS:-} --tag $${CONTAINER_NAME}:${GIT_TAG}"
	elif [ ! -z "$${GIT_BRANCH:-}" ]; then
		export DOCKER_TAGS="$${DOCKER_TAGS:-} --tag $${CONTAINER_NAME}:$${GIT_BRANCH}"
	fi

	docker build --target prod $${DOCKER_TAGS} .

docker-test:
	set -euxo pipefail

	docker run -i --rm \
		--volume="$$(pwd):/workdir" \
		--workdir="/workdir" \
		--user="$$(id -u):$$(id -g)" \
		--ulimit core=0 \
		"$${CONTAINER_NAME}:latest" \
		bash tests/run-cli-tests.sh

docker-push:
	set -euxo pipefail

	# If GIT_TAG is set, then we are about to release to production, so also push the 'latest' tag
	if [ ! -z "$${GIT_TAG:-}" ]; then
		docker push "$${CONTAINER_NAME}:latest"
		docker push "$${CONTAINER_NAME}:${GIT_TAG}"
	elif [ ! -z "$${GIT_BRANCH:-}" ]; then
		docker push "$${CONTAINER_NAME}:$${GIT_BRANCH}"
	fi
