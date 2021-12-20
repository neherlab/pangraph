.PHONY:	all install pangraph environment release documentation clean
.SUFFIXES:
.SECONDARY:

version := 1.7.0

ifeq ($(jc),)
jc := ./vendor/julia-$(version)/bin/julia
endif

jflags := -q --project=.
julia  := julia $(jflags)
srcs   := $(wildcard src/*.jl src/*/*.jl)

datadir   := data/synthetic
testdatum := $(datadir)/test.fa

all:

install: pangraph
	ln -s $$(pwd)/pangraph/bin/pangraph bin/pangraph

environment:
	bin/setup-pangraph

$(datadir):
	mkdir -p $@

$(testdatum): | $(datadir)
	julia $(jflags) -e 'using PanGraph; PanGraph.Simulation.test()'

pangraph: compile.jl trace.jl $(testdatum) $(srcs)
	$(jc) $(jflags) $<

documentation:
	cd docs && julia make.jl

release:
	tar czf pangraph.tar.gz pangraph

clean:
	rm -rf pangraph pangraph.tar.gz

include script/rules.mk
