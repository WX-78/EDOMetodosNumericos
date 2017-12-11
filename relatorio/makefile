compile: command1 command2 command3 command4

command1:
	touch a.bbl a.aux && rm *.bbl *.aux
	gnuplot -e "load 'graficos.gp'"

command2:
	pdflatex --shell-escape relatorio.tex

command3:
	bibtex relatorio

command4:
	pdflatex relatorio.tex
	pdflatex relatorio.tex
command5:
	
	ln -sf imagens/graficoA.tex graficoA.tex
	ln -sf imagens/graficoA.gp  graficoA.gp
	ln -sf imagens/graficoA.eps graficoA.eps
	ln -sf imagens/graficoA.eps graficoB.eps
run:
	evince relatorio.pdf
