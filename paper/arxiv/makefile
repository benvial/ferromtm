

all: clean pdf clean


pdf:
	pdflatex $(file).tex
	bibtex $(file)
	pdflatex $(file).tex
	pdflatex $(file).tex



clean:
		find . \( -name "*.fls" -o -name "*.out" -o -name "*.log" -o -name "*.xdv" \
		-o -name "*.bbl" -o -name "*.blg" -o -name "*.fdb_latexmk" -o -name "*.synctex.gz"\
		-o -name "*.bak" -o -name "*.aux" -o -name "*.xml" -o -name "*.bcf" \
		-o -name "*Notes.bib"  -o -name "*.out"  -o -name "*.fls" \) -exec rm {} \;
