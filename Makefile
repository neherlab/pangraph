DIR := data/synth

.PHONY:	all clean
.SUFFIXES:

all: $(DIR)/guide.json

$(DIR)/seq.fa:
	pangraph generate -d $(DIR) -L 10000 -m 1e-4 -N 50 -T 150

$(DIR)/guide.json: $(DIR)/seq.fa
	pangraph cluster -d $(DIR) $(DIR)/seq.fa

clean:
	@rm -f $(DIR)/*
