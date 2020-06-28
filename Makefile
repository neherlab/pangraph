.PHONY:	all clean
.SUFFIXES:
.SECONDARY:

# set all paths
ROOT_DIR := data/synth
TGT_FILE := targets

DIRS := $(shell sed '/^#.*$$/d' "$(TGT_FILE)")
DIRS := $(addprefix $(ROOT_DIR)/, $(DIRS))

SEQS := $(addsuffix /seq.fa, $(DIRS))
TREE := $(addsuffix /guide.json, $(DIRS))
GRPH := $(addsuffix /graph000.fa, $(DIRS))
# ... etc ...

all: $(GRPH)

%seq.fa:
	@mkdir -p $(@D)
	@$(eval vars=$(subst _, ,$(@D)))
	@$(eval N=$(word 2,$(vars)))
	@$(eval T=$(word 3,$(vars)))
	@$(eval M=$(word 4,$(vars)))
	@$(eval H=$(word 5,$(vars)))
	@$(eval I=$(word 6,$(vars)))
	@$(eval X=$(word 7,$(vars)))
	@echo "generate    "$(@D);\
	pangraph generate -d $(@D) -L 10000 -m $(M) -N $(N) -T $(T) --rate_hgt $(H) --rate_indel $(I) --rate_transpose $(X)

%guide.json: %seq.fa
	@echo "cluster     "$(@D);\
	pangraph cluster -d $(@D) $^

%graph000.fa: %guide.json
	@echo "build	    "$(@D);\
	pangraph build -d $(@D) $^ 2>$(@D)/build.log

clean:
	rm -rf $(ROOT_DIR)/*
