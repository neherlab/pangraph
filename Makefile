DIR := data/synth

.PHONY:	all clean
.SUFFIXES:

all: graph

$(DIR)/seq.fa:
	pangraph generate -d $(DIR)

$(DIR)/tree:
	pangraph cluster -d $(DIR)

clean:
	rm $(DIR)/*
