name=HY_fysiikka_LuKtemplate
bibfile=bibliography.bib

all: $(name).pdf
	
simple: 
	pdflatex $(name).tex
	
clean:
	rm -f *.aux *.blg *.out *.bbl *.log *.nlo *.nls *.toc *.ilg $(name).pdf
	
$(name).pdf: $(name).tex $(bibfile)
	pdflatex $(name).tex
	bibtex $(name).aux
	makeindex $(name).nlo -s nomencl.ist -o $(name).nls
	pdflatex $(name).tex
	pdflatex $(name).tex

