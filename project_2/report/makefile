# LaTeX Makefile
FILE=main
all: $(FILE).pdf
.PHONY: clean
clean:
	\rm *.aux *.blg *.out *.bbl *.log main.pdf
$(FILE).pdf: $(FILE).tex
	pdflatex $(FILE)
	pdflatex $(FILE)
	bibtex $(FILE)
	pdflatex $(FILE)
	pdflatex $(FILE)
