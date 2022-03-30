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
	$(jc) $(jflags) -e 'import Pkg; Pkg.instantiate();'
	$(jc) $(jflags) -e 'import Pkg; Pkg.add(name="Conda"); import Conda; Conda.add("ete3", channel="etetoolkit")' \
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

documentation:
	cd docs && julia --project=./.. make.jl

release:
	tar czf pangraph.tar.gz pangraph

clean:
	rm -rf pangraph pangraph.tar.gz

include script/rules.mk


export CONTAINER_NAME=neherlab/pangraph

SHELL=bash
.ONESHELL:
docker:
	set -euxo pipefail

	# If $RELEASE_VERSION is set, use it as an additional docker tag
	export DOCKER_TAGS="--tag $${CONTAINER_NAME}:latest"
	if [ ! -z "$${RELEASE_VERSION:-}" ]; then
		export DOCKER_TAGS="$${DOCKER_TAGS} --tag $${CONTAINER_NAME}:$${RELEASE_VERSION}"
	fi

	docker build --target prod $${DOCKER_TAGS} .

docker-push:
	set -euxo pipefail
	: "$${RELEASE_VERSION:?The RELEASE_VERSION environment variable is required.}"
	docker push ${CONTAINER_NAME}:${RELEASE_VERSION}
	docker push ${CONTAINER_NAME}:latest
