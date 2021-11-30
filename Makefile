.PHONY:	all pangraph environment documentation clean
.SUFFIXES:
.SECONDARY:

version := 1.6.4

ifeq ($(jc),)
jc := ./vendor/julia-$(version)/bin/julia
endif

jflags := -q --project=.
srcs   := $(wildcard src/*.jl src/*/*.jl)

all: pangraph

environment:
	bin/setup-pangraph

pangraph: compile.jl trace.jl $(srcs)
	$(jc) $(jflags) $<

documentation:
	cd docs && julia make.jl

clean:
	rm -rf pangraph
