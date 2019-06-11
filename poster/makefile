all: clean pdf


pdf:
	xelatex main_poster.tex
	biber main_poster
	xelatex main_poster.tex



clean:
		find . \( -name "*.fls" -o -name "*.out" -o -name "*.log" -o -name "*.xdv" \
		-o -name "*.bbl" -o -name "*.blg" -o -name "*.fdb_latexmk" -o -name "*.synctex.gz"\
		-o -name "*.bak" -o -name "*.aux" -o -name "*.xml" -o -name "*.bcf" \
	 	-o -name "*.toc" -o -name "*.nav" -o -name "*.snm" \
		-o -name "*Notes.bib"  -o -name "*.out"  -o -name "*.fls" \) -exec rm {} \;
