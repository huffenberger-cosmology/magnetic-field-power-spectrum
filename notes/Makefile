magnetic_field_power_spectra.pdf : magnetic_field_power_spectra.tex *.bib
	pdflatex -halt-on-error $<
	bibtex `basename -s .tex $<`
