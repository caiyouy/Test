Latexmk=latexmk

all:report.pdf 

report.pdf:report.tex ./sec/*.tex *.bib
	$(Latexmk) report.tex -time -xelatex

preview:report.pdf
	evince report.pdf&

.PHONY:clean realclean
clean:
	$(Latexmk) -c -time
realclean:
	$(Latexmk) -C -time
