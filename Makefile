bfcsnp:
	Rscript --quiet data-raw/plyranges-ex-snparray.R

bfcbam:
	Rscript --quiet data-raw/plyranges-ex-bam.R

hub:
	Rscript --quiet data-raw/hub.R

pdfbook:
	Rscript --quiet _render.R "bookdown::pdf_book"

gitbook:
	Rscript --quiet _render.R "bookdown::gitbook"

skim:
	open docs/thesis.pdf

preview:
	open docs/index.html

both:
	Rscript --quiet _render.R
