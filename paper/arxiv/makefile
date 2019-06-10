

all: pdf clean


pdf:
	pdflatex $(file).tex
	bibtex $(file)
	pdflatex $(file).tex
	pdflatex $(file).tex



clean:
	rm *.fls*.out *.log *.xdv *.bbl *.blg *.fdb_latexmk *.synctex.gz *.bak *.aux *.xml *.bcf *Notes.bib *.out *.fls
# clean:
# 	find . \( -name "*.fls" -o -name "*.out" -o -name "*.log" -o -name "*.xdv" -o -name "*.bbl" -o -name "*.blg" -o -name "*.fdb_latexmk" -o -name "*.synctex.gz" -o -name "*.bak" -o -name "*.aux" -o -name "*.xml" -o -name "*.bcf" \) -exec rm {} \;
