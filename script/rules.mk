.PHONY: fig1 panX panX2

panx-dir  := data/panx/kleb
fig1-data := $(datadir)/alignment-compare.jld2
benchmark := $(datadir)/benchmark.txt

# comparison to synthetic
$(fig1-data): script/make-sequence.jl script/assay-alignment.jl
	@echo "collecting data for figure 1...";\
	script/make-comparison $^ $@

# benchmark
$(benchmark): script/benchmark
	$^ > $@

# comparison to panX
panx-input-gb := $(wildcard $(panx-dir)/input_GenBank/*.gbk)
panx-input-fa := $(patsubst $(panx-dir)/input_GenBank/%.gbk, $(panx-dir)/fa/%.fa, $(panx-input-gb))

$(panx-dir)/fa:
	mkdir -p $@

$(panx-dir)/fa/%.fa: $(panx-dir)/input_GenBank/%.gbk | $(panx-dir)/fa
	script/gbk-to-fa $* <$< >$@

$(panx-dir)/panx.json: $(panx-input-fa) $(panx-input-gb)
	@echo making $@;\
	script/panx-to-pangraph $(panx-dir)/allclusters_final.tsv $^

$(panx-dir)/pangraph.json: $(panx-input-fa)
	JULIA_NUM_THREADS=8 pangraph build --circular -m 0 -b 50 -l 50 $^ 2>/dev/null 1>$@

panX2: $(panx-dir)/panx.json
panX: $(panx-dir)/pangraph.json

fig1: $(fig1-data)
