# This makes pretty plots of all timing sets.

plotter <- function(data, wrt, col, lty, pch, loc="topleft", cex.axis=1.2, lower=NULL, upper=NULL, ylab="Time (ms)", ...) {
    # Adjusting boundaries. 
    yranges <- range(data$Time)
    if (!is.null(lower)) {
        yranges[1] <- lower
    }
    if (!is.null(upper)) {
        yranges[2] <- upper
    }

    # Creating the empty plot.
    xranges <- range(data[,wrt])
    by.method <- split(data[,c(wrt, "Time")], data$Type, drop=TRUE)
    par(mar=c(5.1, 5.1, 4.1, 2.1))
    plot(1,1,type="n", xlim=xranges, ylim=yranges, log="xy", ylab=ylab, cex.axis=cex.axis, cex.lab=1.4, ..., cex.main=1.4)
 
    # Adding the horizontal tracking lines.
    line.ranges <- range(log2(yranges))
    line.heights <- 2^seq(floor(line.ranges[1]), ceiling(line.ranges[2]), by=1)
    abline(h=line.heights, lwd=0.5, lty=3, col="grey50")

    # Filling in aesthetic details.
    if (missing(col)) { 
        col <- rep("black", length(by.method)) 
        names(col) <- names(by.method)
    }
    if (missing(lty)) { 
        lty <- rep(1, length(by.method)) 
        names(lty) <- names(by.method)
    }
    if (missing(pch)) { 
        pch <- rep(16, length(by.method)) 
        names(pch) <- names(by.method)
    }

    for (meth in names(by.method)) {
        current <- by.method[[meth]]
        points(current[,wrt], current$Time, col=col[meth], pch=pch[meth], cex=1.4)
        lines(current[,wrt], current$Time, col=col[meth], lwd=2, lty=lty[meth])
    }   
    if (!is.na(loc)) {
        legend(loc, col=col, lty=lty[names(col)], pch=pch[names(col)],
               legend=names(col), lwd=2, cex=1.2, bg="white")
    }
    return(invisible(NULL))
}

##############################
# Defining colours.

.beachmat_primary <- "red"
.beachmat_secondary <- "orange"
.beachmat_tertiary <- "gold2"

.Rcpp <- "black"
.Rcpp_arma <- "grey70"
.Rcpp_eigen <- "grey50"

.beachmat_pch_primary <- 16
.beachmat_pch_secondary <- 17
.beachmat_pch_tertiary <- 1

.Rcpp_pch <- 18
.Rcpp_pch_arma <- 15
.Rcpp_pch_eigen <- 4

##############################
# Base plots.

base_cols <- c(beachmat=.beachmat_primary, 
               `beachmat (no copy)`=.beachmat_secondary,
               Rcpp=.Rcpp)
base_pch <- c(beachmat=.beachmat_pch_primary,
              `beachmat (no copy)`=.beachmat_pch_secondary,
              Rcpp=.Rcpp_pch)

incoming <- read.table("../timings_base_col.txt", header=TRUE, sep="\t")
incoming$Ncells <- incoming$Ncells/1e3
pdf("base_col.pdf")
plotter(incoming, "Ncells", col=base_cols, pch=base_pch, xlab=expression("Number of columns ("*10^3*")"), main="Column access")
dev.off()

incoming <- read.table("../timings_base_row.txt", header=TRUE)
incoming$Ngenes <- incoming$Ngenes/1e3
pdf("base_row.pdf")
plotter(incoming, "Ngenes", col=base_cols, pch=base_pch, xlab=expression("Number of rows ("*10^3*")"), main="Row access", loc=NA)
dev.off()

##############################
# Sparse plots.

sparse_cols <- c(beachmat=.beachmat_primary,
                 `beachmat (no copy)`=.beachmat_secondary,
                 `beachmat (ordinary)`=.beachmat_tertiary,
                 RcppArmadillo=.Rcpp_arma,
                 RcppEigen=.Rcpp_eigen)

sparse_pch <- c(beachmat=.beachmat_pch_primary,
                `beachmat (no copy)`=.beachmat_pch_secondary,
                `beachmat (ordinary)`=.beachmat_pch_tertiary,
                RcppArmadillo=.Rcpp_pch_arma,
                RcppEigen=.Rcpp_pch_eigen)

incoming <- read.table("../timings_sparse_col.txt", header=TRUE, sep="\t")
incoming$Ncells <- incoming$Ncells/1e3
incoming$Density <- incoming$Density * 100

pdf("sparse_col_density.pdf")
subincoming <- incoming[incoming$Ncells==1,]
plotter(subincoming, "Density", col=sparse_cols, pch=sparse_pch, xlab="Density (%)", loc="bottomright")
dev.off()

pdf("sparse_col_ncol.pdf")
subincoming <- incoming[incoming$Density==1,]
plotter(subincoming, "Ncells", col=sparse_cols, pch=sparse_pch, xlab=expression("Number of columns ("*10^3*")"), loc=NA)
dev.off()

# By row.

sparse_cols <- c(`beachmat (cached)`=.beachmat_primary,
                 `beachmat (ordinary)`=.beachmat_tertiary,
                 `Rcpp (naive)`=.Rcpp,
                 RcppArmadillo=.Rcpp_arma,
                 RcppEigen=.Rcpp_eigen)

sparse_pch <- c(`beachmat (cached)`=.beachmat_pch_primary,
                `beachmat (ordinary)`=.beachmat_pch_tertiary,
                `Rcpp (naive)`=.Rcpp_pch,
                RcppArmadillo=.Rcpp_pch_arma,
                RcppEigen=.Rcpp_pch_eigen)

incoming <- read.table("../timings_sparse_row.txt", header=TRUE, stringsAsFactors=FALSE, sep="\t")
incoming$Ngenes <- incoming$Ngenes/1e3
incoming$Density <- incoming$Density * 100

pdf("sparse_row_density.pdf")
subincoming <- incoming[incoming$Ngenes==10,]
plotter(subincoming, "Density", col=sparse_cols, pch=sparse_pch, xlab="Density (%)", upper=max(subincoming$Time)*1.3)
dev.off()

pdf("sparse_row_nrow.pdf")
subincoming <- incoming[incoming$Density==1,]
plotter(subincoming, "Ngenes", col=sparse_cols, pch=sparse_pch, xlab=expression("Number of rows ("*10^3*")"), loc=NA)
dev.off()

# By non-consecutive rows.

sparse_cols <- c(`beachmat (cached)`=.beachmat_primary,
                 `beachmat (ordinary)`=.beachmat_tertiary,
                 `Rcpp (naive)`=.Rcpp)

sparse_pch <- c(`beachmat (cached)`=.beachmat_pch_primary,
                `beachmat (ordinary)`=.beachmat_pch_tertiary,
                `Rcpp (naive)`=.Rcpp_pch)

incoming <- read.table("../timings_sparse_row_ordered.txt", header=TRUE, stringsAsFactors=FALSE, sep="\t")
incoming$Ngenes <- incoming$Ngenes/1e3
pdf("sparse_row_ordered.pdf")
plotter(incoming, "Ngenes", col=sparse_cols, pch=sparse_pch, xlab=expression("Number of rows ("*10^3*")"), main="Ordered")
dev.off()

incoming <- read.table("../timings_sparse_row_random.txt", header=TRUE, stringsAsFactors=FALSE, sep="\t")
incoming$Ngenes <- incoming$Ngenes/1e3
pdf("sparse_row_random.pdf")
plotter(incoming, "Ngenes", col=sparse_cols, pch=sparse_pch, xlab=expression("Number of rows ("*10^3*")"), main="Random", loc=NA)
dev.off()

##############################
# Making HDF5 plots.

hdf5_cols <- c(`HDF5 (rectangle)`=.beachmat_primary,
               `HDF5 (column)`=.beachmat_secondary,
               `ordinary`=.beachmat_tertiary)

hdf5_pch <- c(`HDF5 (rectangle)`=.beachmat_pch_primary,
              `HDF5 (column)`=.beachmat_pch_secondary,
              `ordinary`=.beachmat_pch_tertiary)

incoming <- read.table("../timings_hdf5_col.txt", header=TRUE, sep="\t")
incoming$Ncells <- incoming$Ncells/1e3

pdf("HDF5_col_ncol.pdf")
plotter(incoming, "Ncells", hdf5_cols, pch=hdf5_pch, xlab=expression("Number of columns ("*10^3*")"), 
        main="Column access", cex.axis=1, lower=16, upper=16000)
dev.off()

incoming <- read.table("../timings_hdf5_row.txt", header=TRUE, sep="\t")
incoming$Ngenes <- incoming$Ngenes/1e3

names(hdf5_cols) <- sub("column", "row", names(hdf5_cols))
names(hdf5_pch) <- sub("column", "row", names(hdf5_pch))

pdf("HDF5_row_nrow.pdf")
plotter(incoming, "Ngenes", hdf5_cols, pch=hdf5_pch, xlab=expression("Number of rows ("*10^3*")"), 
        main="Row access", cex.axis=1, lower=16, upper=16000)
dev.off()

##############################
# Layout type; this involves a separate colour scheme.

layout_cols <- c(Contiguous="darkolivegreen",
                 `Column chunks`="forestgreen",
                 `Row chunks`="limegreen",
                 `Rectangular chunks`="yellowgreen")
layout_pch <- c(Contiguous=17,
                `Column chunks`=15,
                `Row chunks`=18,
                `Rectangular chunks`=16)

# By column:

incoming <- read.table("../timings_hdf5_col_layout.txt", header=TRUE, sep="\t")
incoming <- incoming[-grep("uncompressed", incoming$Type),]

pdf("HDF5_col_layout.pdf")
plotter(incoming, "Ncells", layout_cols, pch=layout_pch,
    xlab="Number of columns", main="Column access", cex.axis=1, upper=50000)
dev.off()

# By row:

incoming <- read.table("../timings_hdf5_row_layout.txt", header=TRUE, sep="\t")
incoming <- incoming[-grep("uncompressed", incoming$Type),]
incoming$Time <- incoming$Time/1e3

pdf("HDF5_row_layout.pdf")
plotter(incoming, "Ngenes", layout_cols, pch=layout_pch,
    xlab="Number of rows", main="Row access", cex.axis=1, ylab="Time (s)", yaxt="n", loc=NA)
ticks <- c(10^(-2:2))
axis(2, at=ticks, ticks)
dev.off()

# By column (random):

incoming <- read.table("../timings_hdf5_col_layout_random.txt", header=TRUE, sep="\t")

pdf("HDF5_col_layout_random.pdf")
plotter(incoming, "Ncells", layout_cols, pch=layout_pch,
    xlab="Number of columns", main="Column access", cex.axis=1, upper=max(incoming$Time)*4)
dev.off()

# By row (random):

incoming <- read.table("../timings_hdf5_row_layout_random.txt", header=TRUE, sep="\t")

pdf("HDF5_row_layout_random.pdf")
plotter(incoming, "Ngenes", layout_cols, pch=layout_pch,
    xlab="Number of rows", main="Row access", cex.axis=1, loc=NA)
dev.off()

##############################
# Rechunking

incoming <- read.table("../chunking/timings_rechunk_col.txt", header=TRUE, sep="\t")
incoming$Ncells <- incoming$Ncells/1e3

cols <- c(`Column to row`="darkgreen",
          `Row to column`="limegreen")

pchs <- c(`Column to row`=16,
          `Row to column`=17)

pdf("HDF5_col_rechunk.pdf")
plotter(incoming, "Ncells", col=cols, pch=pchs,
        xlab=expression("Number of columns ("*10^3*")"), cex.axis=1)
dev.off()

incoming <- read.table("../chunking/timings_rechunk_row.txt", header=TRUE, sep="\t")
incoming$Ngenes <- incoming$Ngenes/1e3

pdf("HDF5_row_rechunk.pdf")
plotter(incoming, "Ngenes", col=cols, pch=pchs,
        xlab=expression("Number of rows ("*10^3*")"), cex.axis=1, loc=NA)
dev.off()

##############################
# Matrix multiplication. Again, this involves a different suite of colors. 

mult_cols <- c(`ordinary (beachmat)`="darkblue",
               `ordinary (R)`="darkblue",
               `sparse (beachmat)`="blue",
               `sparse (beachmat II)`="blue",
               `sparse (R)`="blue",
               `HDF5/ordinary (beachmat)`="dodgerblue",
               `HDF5/ordinary (R)`="dodgerblue",
               `HDF5/HDF5 (beachmat)`="dodgerblue")

mult_pch <- c(`ordinary (beachmat)`=16,
              `ordinary (R)`=17,
              `sparse (beachmat)`=18,
              `sparse (beachmat II)`=15,
              `sparse (R)`=4,
              `HDF5/ordinary (beachmat)`=1,
              `HDF5/HDF5 (beachmat)`=2,
              `HDF5/ordinary (R)`=3)

mult_lty <- c(`ordinary (beachmat)`=1,
              `ordinary (R)`=2,
              `HDF5/ordinary (beachmat)`=1,
              `HDF5/HDF5 (beachmat)`=1,
              `HDF5/ordinary (R)`=2,
              `sparse (beachmat)`=1,
              `sparse (beachmat II)`=1,
              `sparse (R)`=2)

incoming <- read.table("../timings_mat_mult.txt", header=TRUE, sep="\t")

pdf("mat_mult.pdf", width=11, height=8)
layout(cbind(1, 2), width=c(5, 2))
par(mar=c(5.1, 4.1, 2.1, 0.1))
plotter(incoming, "N", mult_cols, pch=mult_pch, lty=mult_lty,
        xlab=expression("Order ("*10^3*")"), cex.axis=1, loc=NA)
par(mar=c(5.1, 0.1, 2.1, 0.1))
plot.new()
legend("left", col=mult_cols, lty=mult_lty[names(mult_cols)], pch=mult_pch[names(mult_cols)],
       legend=names(mult_cols), lwd=2, cex=1.2, bg="white")
dev.off()
