.PHONY:	all install pangraph environment release documentation clean
.SUFFIXES:
.SECONDARY:

version := 1.6.4

ifeq ($(jc),)
jc := ./vendor/julia-$(version)/bin/julia
endif

jflags := -q --project=.
srcs   := $(wildcard src/*.jl src/*/*.jl)

all: install

install: pangraph
	ln -s $(pwd)/pangraph/bin/pangraph bin/$<

environment:
	bin/setup-pangraph

pangraph: compile.jl trace.jl $(srcs)
	$(jc) $(jflags) $<

documentation:
	cd docs && julia make.jl

release:
	tar czf pangraph.tar.gz pangraph

clean:
	rm -rf pangraph pangraph.tar.gz
