.PHONY: fig1 panX

panx-dir  := data/panx/kleb
fig1-data := $(datadir)/alignment-compare.jld2

# comparison to synthetic
$(fig1-data): script/make-sequence.jl script/assay-alignment.jl
	@echo "collecting data for figure 1...";\
	script/make-comparison $^ $@

# comparison to panX
panx-input := $(patsubst $(panx-dir)/input_GenBank/%.gbk, $(panx-dir)/fa/%.fa, $(wildcard $(panx-dir)/input_GenBank/*.gbk))

$(panx-dir)/fa:
	mkdir -p $@

$(panx-dir)/fa/%.fa: $(panx-dir)/input_GenBank/%.gbk | $(panx-dir)/fa
	script/gbk-to-fa $* <$< >$@

$(panx-dir)/pangraph.json: $(panx-input)
	pangraph build --circular -m 0 -b 50 -l 50 $^ 2>/dev/null 1>$@

panX: $(panx-dir)/pangraph.json

fig1: $(fig1-data)
