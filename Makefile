.PHONY:	all install pangraph environment release documentation clean
.SUFFIXES:
.SECONDARY:

version := 1.7.2

ifeq ($(jc),)
jc := ./vendor/julia-$(version)/bin/julia
endif

jflags := -q --project=.
julia  := julia $(jflags)
srcs   := $(wildcard src/*.jl src/*/*.jl)

datadir   := data/synthetic
testdatum := $(datadir)/test.fa

all: pangraph install

install: pangraph/bin/pangraph
	ln -s $$(pwd)/$< bin/pangraph

environment:
	bin/setup-pangraph

pangraph: pangraph/bin/pangraph

$(datadir):
	mkdir -p $@

$(testdatum): | $(jc) $(datadir)
	$(jc) $(jflags) -e 'import Pkg; Pkg.instantiate(); Pkg.build()'
	$(jc) $(jflags) -e 'using PanGraph; PanGraph.Simulation.test()'

# TODO: look for ARM vs x86
$(jc):
#ifeq ($(shell uname -s),Linux)
#	cd vendor && \
#	curl -L https://julialang-s3.julialang.org/bin/linux/x64/$(basename $(version))/julia-$(version)-linux-x86_64.tar.gz -o julia-$(version)-linux-x86_64.tar.gz && \
#	tar xzf julia-$(version)-linux-x86_64.tar.gz
#else
#ifeq ($(shell uname -s),Darwin)
#	cd vendor && \
#	curl -L https://julialang-s3.julialang.org/bin/mac/x64/$(basename $(version))/julia-$(version)-mac64.tar.gz -o julia-$(version)-mac64.tar.gz && \
#	tar xzf julia-$(version)-mac64.tar.gz
#else
#	$(error unsupported host system)
#endif
#endif

pangraph/bin/pangraph: compile.jl trace.jl $(srcs) $(testdatum) $(jc)
	$(jc) $(jflags) $<

documentation:
	cd docs && julia --project=./.. make.jl

release:
	tar czf pangraph.tar.gz pangraph

clean:
	rm -rf pangraph pangraph.tar.gz

include script/rules.mk

SHELL=bash
.ONESHELL:
docker:
	set -euxo pipefail

	export CONTAINER_NAME=neherlab/pangraph

	docker build \
	--target prod \
	--tag $${CONTAINER_NAME} \
	--build-arg UID=$(shell id -u) \
	--build-arg GID=$(shell id -g) \
	.

	docker run -it --rm \
	--name=$${CONTAINER_NAME}-$(shell date +%s) \
	--init \
	--user=$(shell id -u):$(shell id -g) \
	--volume="$(shell pwd):/workdir" \
	--workdir="/workdir" \
	$${CONTAINER_NAME}
