all: template.pdf clean

template.pdf:template.tex 
	pdflatex template.tex

.PHONY:clean realclean
clean:
	-rm *.aux *.dvi *.log

realclean:clean
	-rm *.pdf
