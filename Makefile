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
