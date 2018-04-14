all: Figure1.pdf Figure2.pdf Figure3.pdf Figure4.pdf Figure5.pdf S1_Fig.pdf S2_Fig.pdf S3_Fig.pdf S4_Fig.pdf S5_Fig.pdf S6_Fig.pdf S7_Fig.pdf S8_Fig.pdf

Figure1.pdf: pics/overview.svg
	inkscape -D -A $@ $<

Figure2.pdf: pics/figure_sparse_times.tex ../simulations/timings/pics/sparse_row_density.pdf ../simulations/timings/pics/sparse_row_nrow.pdf
	cd pics; pdflatex `basename $<`
	pdfcrop pics/figure_sparse_times.pdf $@

Figure3.pdf: pics/figure_hdf5_times.tex ../simulations/timings/pics/HDF5_col_layout.pdf ../simulations/timings/pics/HDF5_row_layout.pdf
	cd pics; pdflatex `basename $<`
	pdfcrop pics/figure_hdf5_times.pdf $@

Figure4.pdf: pics/figure_zeisel.tex ../real/zeisel/pics/zeisel_col.pdf ../real/zeisel/pics/zeisel_row.pdf ../real/zeisel/pics/zeisel_detection.pdf
	cd pics; pdflatex `basename $<`
	pdfcrop pics/figure_zeisel.pdf $@

Figure5.pdf: pics/figure_10x.tex ../real/10X/pics/cycle.pdf ../real/10X/pics/sizefacs.pdf ../real/10X/pics/hvg.pdf ../real/10X/pics/pca.pdf
	cd pics; pdflatex `basename $<`
	pdfcrop pics/figure_10x.pdf $@

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
   
