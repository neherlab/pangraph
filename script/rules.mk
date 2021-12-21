.PHONY: fig1 panx

panx-dir  := data/panx/kleb
accuracy  := $(datadir)/alignment-compare.jld2
benchmark := $(datadir)/benchmark.txt

# benchmarking
$(benchmark): script/benchmark
	$^ > $@

# accuracy relative to simulated
$(accuracy): script/make-sequence.jl script/make-accuracy.jl
	@echo "collecting data for figure 1...";\
	script/make-accuracy $^ $@

# accuracy relative to panX gene boundaries
panx-input-gb := $(wildcard $(panx-dir)/input_GenBank/*.gbk)
panx-input-fa := $(patsubst $(panx-dir)/input_GenBank/%.gbk, $(panx-dir)/fa/%.fa, $(panx-input-gb))

$(panx-dir)/fa:
	mkdir -p $@

$(panx-dir)/fa/%.fa: $(panx-dir)/input_GenBank/%.gbk | $(panx-dir)/fa
	@echo MAKE	$@;\
	script/gbk-to-fa $* <$< >$@

$(panx-dir)/panx.json: $(panx-input-fa) $(panx-input-gb)
	@echo MAKE	$@;\
	script/panx-to-pangraph $(panx-dir)/vis/geneCluster.json $@ $^

$(panx-dir)/pangraph.json: $(panx-input-fa)
	@echo MAKE	$@;\
	JULIA_NUM_THREADS=8 pangraph build --circular -m 0 -b 0 -l 50 $^ 1>$@

panx: $(panx-dir)/pangraph.json $(panx-dir)/panx.json

fig1: $(accuracy)
