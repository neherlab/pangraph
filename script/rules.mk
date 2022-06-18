.PHONY: fig1 panx

accuracy  := $(datadir)/accuracy
benchmark := $(datadir)/benchmark.txt

# ------------------------------------------------------------------------
#  synthetic data

# benchmarking
$(benchmark): script/benchmark
	@echo "MAKE $@";\
	$^ > $@

# accuracy relative to simulated
$(accuracy)-%.jld2: script/make-sequence.jl script/make-accuracy.jl
	@echo "MAKE $@";\
	script/make-accuracy $* $^ $@

figs/benchmark.png: $(benchmark)
	julia --project=script script/plot-benchmark.jl $< $@

figs/paper-accuracy-%.png: $(accuracy)-%.jld2
	julia --project=script script/plot-accuracy.jl $< figs

# ------------------------------------------------------------------------
#  comparison to panx for multiple species

# accuracy relative to panX gene boundaries
define PANX
$(eval input-gb := $(wildcard script/panx_data/$(1)/input_GenBank/*.gbk))
$(eval input-fa := $(patsubst script/panx_data/$(1)/input_GenBank/%.gbk, script/panx_data/$(1)/fa/%.fa, $(input-gb)))
script/panx_data/$(1)/fa:
	mkdir -p $$@

script/panx_data/$(1)/fa/%.fa: script/panx_data/$(1)/input_GenBank/%.gbk | script/panx_data/$(1)/fa
	@echo "MAKE	$$@";\
	script/gbk-to-fa $$* <$$< >$$@

script/panx_data/$(1)/panx.json: $(input-fa) $(input-gb)
	@echo "MAKE	$$@";\
	script/panx-to-pangraph script/panx_data/$(1)/vis/geneCluster.json $$@ $$^

script/panx_data/$(1)/pangraph.json: $(input-fa)
	@echo "MAKE	$$@";\
	JULIA_NUM_THREADS=8 pangraph build --circular --upper-case --max-self-map 50 $$^ 1>$$@

$(eval panx-targets += script/panx_data/$(1)/pangraph.json script/panx_data/$(1)/panx.json)
endef

panx-targets =
$(eval $(call PANX,klebsiella_pneumoniae))
$(eval $(call PANX,helicobacter_pylori))
$(eval $(call PANX,escherichia_coli))
$(eval $(call PANX,mycobacterium_tuberculosis))
$(eval $(call PANX,prochlorococcus_marinus))

figs/panx-compare.png: $(panx-targets)
	julia --project=script script/plot-panx-compare.jl $^
