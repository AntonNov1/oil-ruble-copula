MAIN = main
LATEX = pdflatex
BIBTEX = bibtex
FLAGS = -interaction=nonstopmode -halt-on-error

.PHONY: all clean view fast

all:
	$(LATEX) $(FLAGS) $(MAIN)
	$(BIBTEX) $(MAIN)
	$(LATEX) $(FLAGS) $(MAIN)
	$(LATEX) $(FLAGS) $(MAIN)

# Быстрая компиляция без библиографии (для правки текста)
fast:
	$(LATEX) $(FLAGS) $(MAIN)

view:
	open $(MAIN).pdf

clean:
	rm -f *.aux *.log *.bbl *.blg *.toc *.out *.lof *.lot *.fls *.fdb_latexmk
