all: main.pdf \
	Figure1.eps Figure2.eps Figure3.eps Figure4.eps Figure5.eps \
	S1_Fig.pdf S2_Fig.pdf S3_Fig.pdf S4_Fig.pdf S5_Fig.pdf S6_Fig.pdf S7_Fig.pdf S8_Fig.pdf \
	S1_Text.pdf S2_Text.pdf S3_Text.pdf S4_Text.pdf

Figure1.eps: pics/overview.svg
	inkscape -D -E $@ $<

Figure2.eps: pics/figure_sparse_times.tex ../simulations/timings/pics/sparse_row_density.pdf ../simulations/timings/pics/sparse_row_nrow.pdf
	cd pics; pdflatex `basename $<`
	pdfcrop `echo $< | sed "s/tex/pdf/"` tmp.pdf
	pdftops -eps tmp.pdf $@
	rm tmp.pdf

Figure3.eps: pics/figure_hdf5_times.tex ../simulations/timings/pics/HDF5_col_layout.pdf ../simulations/timings/pics/HDF5_row_layout.pdf
	cd pics; pdflatex `basename $<`
	pdfcrop `echo $< | sed "s/tex/pdf/"` tmp.pdf
	pdftops -eps tmp.pdf $@
	rm tmp.pdf

Figure4.eps: pics/figure_zeisel.tex ../real/zeisel/pics/zeisel_col.pdf ../real/zeisel/pics/zeisel_row.pdf ../real/zeisel/pics/zeisel_detection.pdf
	cd pics; pdflatex `basename $<`
	pdfcrop `echo $< | sed "s/tex/pdf/"` tmp.pdf
	pdftops -eps tmp.pdf $@
	rm tmp.pdf

Figure5.eps: pics/figure_10x.tex ../real/10X/pics/cycle.pdf ../real/10X/pics/sizefacs.pdf ../real/10X/pics/hvg.pdf ../real/10X/pics/pca.pdf
	cd pics; pdflatex `basename $<`
	pdfcrop `echo $< | sed "s/tex/pdf/"` tmp.pdf
	pdftops -eps tmp.pdf $@
	rm tmp.pdf

S1_Fig.pdf: pics/supp_simple.tex
	cd pics; pdflatex `basename $<`
	pdfcrop `echo $< | sed 's/tex/pdf/'` $@
 
S2_Fig.pdf: pics/sparse.svg
	inkscape -D -A $@ $<

S3_Fig.pdf: pics/supp_sparse_col.tex
	cd pics; pdflatex `basename $<`
	pdfcrop `echo $< | sed 's/tex/pdf/'` $@

S4_Fig.pdf: pics/supp_sparse_row_random.tex
	cd pics; pdflatex `basename $<`
	pdfcrop `echo $< | sed 's/tex/pdf/'` $@

S5_Fig.pdf: pics/supp_hdf5_consecutive.tex
	cd pics; pdflatex `basename $<`
	pdfcrop `echo $< | sed 's/tex/pdf/'` $@

S6_Fig.pdf: pics/supp_hdf5_random.tex
	cd pics; pdflatex `basename $<`
	pdfcrop `echo $< | sed 's/tex/pdf/'` $@

S7_Fig.pdf: pics/supp_hdf5_rechunk.tex
	cd pics; pdflatex `basename $<`
	pdfcrop `echo $< | sed 's/tex/pdf/'` $@

S8_Fig.pdf: ../simulations/timings/pics/mat_mult.pdf
	pdfcrop $< $@

S1_Text.pdf: S1_text.tex
	pdflatex $<
	bibtex `echo $< | sed "s/\\.tex//"`
	pdflatex $<
	pdflatex $<
	mv `echo $< | sed "s/\\.tex/.pdf/"` $@

S2_Text.pdf: S2_text.tex
	pdflatex $<
	bibtex `echo $< | sed "s/\\.tex//"`
	pdflatex $<
	pdflatex $<
	mv `echo $< | sed "s/\\.tex/.pdf/"` $@

S3_Text.pdf: S3_text.tex
	pdflatex $<
	mv `echo $< | sed "s/\\.tex/.pdf/"` $@

S4_Text.pdf: S4_text.tex
	pdflatex $<
	bibtex `echo $< | sed "s/\\.tex//"`
	pdflatex $<
	pdflatex $<
	mv `echo $< | sed "s/\\.tex/.pdf/"` $@

main.pdf: description.tex
	pdflatex $<
	bibtex `echo $< | sed "s/\\.tex//"`
	pdflatex $<
	pdflatex $<
	echo "X <- as.list(readLines('description.tex')); bib.loc <- grep('bibliography.ref.', X); supp.loc <- grep('input.supp_refs.tex.', X); X[[bib.loc]] <- readLines('description.bbl'); X[[supp.loc]] <- readLines('supp_refs.tex'); writeLines(unlist(X), con='main.tex')" | R --slave
	pdflatex main.tex
	pdflatex main.tex
