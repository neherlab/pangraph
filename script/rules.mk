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
$(eval input-gb := $(wildcard data/panx/$(1)/input_GenBank/*.gbk))
$(eval input-fa := $(patsubst data/panx/$(1)/input_GenBank/%.gbk, data/panx/$(1)/fa/%.fa, $(input-gb)))
data/panx/$(1)/fa:
	mkdir -p $$@

data/panx/$(1)/fa/%.fa: data/panx/$(1)/input_GenBank/%.gbk | data/panx/$(1)/fa
	@echo "MAKE	$$@";\
	script/gbk-to-fa $$* <$$< >$$@

data/panx/$(1)/panx.json: $(input-fa) $(input-gb)
	@echo "MAKE	$$@";\
	script/panx-to-pangraph data/panx/$(1)/vis/geneCluster.json $$@ $$^

data/panx/$(1)/pangraph.json: $(input-fa)
	@echo "MAKE	$$@";\
	JULIA_NUM_THREADS=8 pangraph build --circular --upper-case --max-self-map 50 $$^ 1>$$@

$(eval panx-targets += data/panx/$(1)/pangraph.json data/panx/$(1)/panx.json)
endef

panx-targets =
$(eval $(call PANX,kleb))
$(eval $(call PANX,helio))
$(eval $(call PANX,ecoli))
$(eval $(call PANX,myco))
$(eval $(call PANX,proc))

figs/panx-compare.png: $(panx-targets)
	julia --project=script script/plot-panx-compare.jl $^
