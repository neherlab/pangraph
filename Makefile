.PHONY:	all app clean
.SUFFIXES:
.SECONDARY:

ifeq ($(JC),)
JC := ./vendor/julia-1.6.0/bin/julia
endif

JFLAGS := -q --project=.

all: pangraph

pangraph:
	$(JC) $(JFLAGS) compile.jl

clean:
	rm -rf pangraph

