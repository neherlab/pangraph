.PHONY: fig1 panx

panx-dir  := data/panx/kleb
accuracy  := $(datadir)/accuracy
benchmark := $(datadir)/benchmark.txt

# benchmarking
$(benchmark): script/benchmark
	@echo "MAKE $@";\
	$^ > $@

# accuracy relative to simulated
$(accuracy)-%.jld2: script/make-sequence.jl script/make-accuracy.jl
	@echo "MAKE $@";\
	script/make-accuracy $* $^ $@

# accuracy relative to panX gene boundaries
panx-input-gb := $(wildcard $(panx-dir)/input_GenBank/*.gbk)
panx-input-fa := $(patsubst $(panx-dir)/input_GenBank/%.gbk, $(panx-dir)/fa/%.fa, $(panx-input-gb))

$(panx-dir)/fa:
	mkdir -p $@

$(panx-dir)/fa/%.fa: $(panx-dir)/input_GenBank/%.gbk | $(panx-dir)/fa
	@echo "MAKE	$@";\
	script/gbk-to-fa $* <$< >$@

$(panx-dir)/panx.json: $(panx-input-fa) $(panx-input-gb)
	@echo "MAKE	$@";\
	script/panx-to-pangraph $(panx-dir)/vis/geneCluster.json $@ $^

$(panx-dir)/pangraph.json: $(panx-input-fa)
	@echo "MAKE	$@";\
	JULIA_NUM_THREADS=8 pangraph build --circular --upper-case --max-self-map 50 $^ 1>$@

panx: $(panx-dir)/pangraph.json $(panx-dir)/panx.json

fig1: $(accuracy)-10.jld2 $(accuracy)-20.jld2 $(benchmark)
