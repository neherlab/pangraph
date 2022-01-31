.PHONY:	all install pangraph environment release documentation clean
.SUFFIXES:
.SECONDARY:

version := 1.7.1

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

$(testdatum): | $(datadir)
	julia $(jflags) -e 'using PanGraph; PanGraph.Simulation.test()'

$(jc):
	cd vendor && \
	curl -L https://julialang-s3.julialang.org/bin/linux/x64/$(basename $(version))/julia-$(version)-linux-x86_64.tar.gz -o julia-$(version)-linux-x86_64.tar.gz && \
	tar xzf julia-$(version)-linux-x86_64.tar.gz

pangraph/bin/pangraph: compile.jl trace.jl $(testdatum) $(srcs) $(jc)
	$(jc) $(jflags) $<

documentation:
	cd docs && julia make.jl

release:
	tar czf pangraph.tar.gz pangraph

clean:
	rm -rf pangraph pangraph.tar.gz

include script/rules.mk
