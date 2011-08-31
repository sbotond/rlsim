# Makefile for building rlsim documentation.
name=rlsim_manual

$(name).pdf: ref.bib $(name).tex *.tex pix/*
	pdflatex $(name); bibtex $(name); pdflatex $(name); pdflatex $(name)

clean:
	rm -f $(name).log $(name).aux $(name).pdf $(name).out $(name).bbl $(name).blg *.toc *.bak

bib:
	kbibtex ref.bib

v: $(name).pdf
	evince $(name).pdf

com:
	git commit -a

