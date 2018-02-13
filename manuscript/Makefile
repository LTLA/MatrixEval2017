all: Figure2.pdf  Figure3.pdf

Figure2.pdf: pics/figure_sparse_times.tex 
	cd pics; pdflatex `basename $<`
	pdfcrop pics/figure_sparse_times.pdf $@

Figure3.pdf: pics/figure_hdf5_times.tex 
	cd pics; pdflatex `basename $<`
	pdfcrop pics/figure_hdf5_times.pdf $@

