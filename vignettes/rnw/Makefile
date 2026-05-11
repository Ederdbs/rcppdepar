
sources :=              $(wildcard *.Rnw)
vignettes :=            $(sources:.Rnw=.pdf)

all:                    ${vignettes}

%.pdf:                  %.Rnw
			Rscript -e 'utils::Sweave("$^")'
			pdflatex RcppDE.tex
			bibtex RcppDE.aux
			pdflatex RcppDE.tex
			pdflatex RcppDE.tex
			Rscript -e 'tools::compactPDF("$@", gs_quality="ebook")'
			mv -v $@ ..

clean:
			rm -rf *.tex *.log *.aux *.bbl *.blg *.out *~ auto/ figures/
