.PHONY: fig1

fig1-data := $(datadir)/alignment-compare.jld2

$(fig1-data): script/make-sequence.jl script/assay-alignment.jl
	@echo "collecting data for figure 1...";\
	script/make-comparison $^ $@

fig1: $(fig1-data)
