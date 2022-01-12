TEX   := latexmk
FLAGS := -pdf -quiet -bibtex -use-make -pdflatex="lualatex -interaction=nonstopmode"

all: paper/draft.pdf

%.pdf: %.tex
	$(TEX) $(FLAGS) -outdir=$(@D) $<

clean:
	$(TEX) -quiet -cd -C paper/draft.tex
