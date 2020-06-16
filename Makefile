.PHONY:	all clean
.SUFFIXES:

# set all paths
ROOT_DIR := data/synth
TGT_FILE := targets

DIRS != cat "$(TGT_FILE)"
DIRS := $(addprefix $(ROOT_DIR)/, $(DIRS))

SEQS := $(addsuffix /seq.fa, $(DIRS))
TREE := $(addsuffix /guide.json, $(DIRS))
# ... etc ...

all: $(TREE)

%seq.fa:
	@mkdir -p $(@D)
	# extract variables
	@$(eval vars=$(subst _, ,$(@D)))
	@$(eval N=$(word 2,$(vars)))
	@$(eval T=$(word 3,$(vars)))
	@$(eval H=$(word 4,$(vars)))
	@$(eval I=$(word 5,$(vars)))
	@$(eval X=$(word 6,$(vars)))
	pangraph generate -d $(@D) -L 10000 -m 1e-4 -N $(N) -T $(T) --rate_hgt $(H) --rate_indel $(I) --rate_transpose $(X)

%guide.json: %seq.fa
	pangraph cluster -d $(@D) $^

clean:
	@rm -rf $(ROOT_DIR)/*
