# This makes pretty plots of all timing sets.

plotter <- function(data, wrt, col, lty, pch, loc="topleft", cex.axis=1.2, upper=NULL, ylab="Time (ms)", ...) {
    yranges <- range(data$Time)
    if (!is.null(upper)) {
        yranges[2] <- upper        
    }
    xranges <- range(data[,wrt])
    by.method <- split(data[,c(wrt, "Time")], data$Type, drop=TRUE)
    par(mar=c(5.1, 5.1, 4.1, 2.1))
    plot(1,1,type="n", xlim=xranges, ylim=yranges, log="xy", ylab=ylab, cex.axis=cex.axis, cex.lab=1.4, ..., cex.main=1.4)
 
    line.ranges <- range(log2(yranges))
    line.heights <- 2^seq(floor(line.ranges[1]), ceiling(line.ranges[2]), by=1)
    abline(h=line.heights, lwd=0.5, lty=3, col="grey50")

    if (missing(col)) { col <- rep("black", length(by.method)) }
    if (missing(lty)) { lty <- rep(1, length(by.method)) }
    if (missing(pch)) { pch <- rep(16, length(by.method)) }

    for (i in seq_along(by.method)) {
        current <- by.method[[i]]
        points(current[,wrt], current$Time, col=col[i], pch=pch[i], cex=1.4)
        lines(current[,wrt], current$Time, col=col[i], lwd=2, lty=lty[i])
    }   
    if (!is.na(loc)) {
        legend(loc, col=col, lty=lty, legend=names(by.method), pch=pch, lwd=2, cex=1.2, bg="white")
    }
    return(invisible(NULL))
}

##############################
# Base plots.

incoming <- read.table("../timings_base_col.txt", header=TRUE)
incoming$Ncells <- incoming$Ncells/1e3
pdf("base_col.pdf")
plotter(incoming, "Ncells", c("black", "grey70"), pch=c(16, 17), xlab=expression("Number of columns ("*10^3*")"), main="Column access times")
dev.off()

incoming <- read.table("../timings_base_row.txt", header=TRUE)
incoming$Ngenes <- incoming$Ngenes/1e3
pdf("base_row.pdf")
plotter(incoming, "Ngenes", c("black", "grey70"), pch=c(16, 17), xlab=expression("Number of rows ("*10^3*")"), main="Row access times", loc=NA)
dev.off()

##############################
# Sparse plots.

incoming <- read.table("../timings_sparse_col.txt", header=TRUE)
incoming$Ncells <- incoming$Ncells/1e3
incoming$Density <- incoming$Density * 100

pdf("sparse_col_density.pdf")
subincoming <- incoming[incoming$Ncells==1,]
plotter(subincoming, "Density", c("grey70", "black", "red"), pch=c(18, 16, 17), xlab="Density (%)", loc="bottomright")
dev.off()

pdf("sparse_col_ncol.pdf")
subincoming <- incoming[incoming$Density==1,]
plotter(subincoming, "Ncells", c("grey70", "black", "red"), pch=c(18, 16, 17), xlab=expression("Number of columns ("*10^3*")"), loc=NA)
dev.off()

# By row.

incoming <- read.table("../timings_sparse_row.txt", header=TRUE)
incoming$Ngenes <- incoming$Ngenes/1e3
incoming$Density <- incoming$Density * 100
incoming <- incoming[incoming$Type!="arma",]

pdf("sparse_row_density.pdf")
subincoming <- incoming[incoming$Ngenes==10,]
plotter(subincoming, "Density", c("black", "red"), pch=c(16, 17), xlab="Density (%)")
dev.off()

pdf("sparse_row_nrow.pdf")
subincoming <- incoming[incoming$Density==1,]
plotter(subincoming, "Ngenes", c("black", "red"), pch=c(16, 17), xlab=expression("Number of rows ("*10^3*")"), loc=NA)
dev.off()

# By access pattern.

incoming <- read.table("../timings_sparse_row_naive.txt", header=TRUE)
incoming$Ngenes <- incoming$Ngenes/1e3
incoming$Density <- incoming$Density * 100

pdf("sparse_row_naive_density.pdf")
subincoming <- incoming[incoming$Ngenes==10,]
plotter(subincoming, "Density", c("red", "blue"), pch=c(17, 15), lty=c(2, 2), xlab="Density (%)")
dev.off()

pdf("sparse_row_naive_nrow.pdf")
subincoming <- incoming[incoming$Density==1,]
plotter(subincoming, "Ngenes", c("red", "blue"), pch=c(17, 15), lty=c(2,2), xlab=expression("Number of rows ("*10^3*")"), loc=NA)
dev.off()

##############################
# Making HDF5 plots.

incoming <- read.table("../timings_hdf5_col.txt", header=TRUE)
incoming$Ncells <- incoming$Ncells/1e3

pdf("HDF5_col_ncol.pdf")
plotter(incoming, "Ncells", c("red", "black"), pch=c(18, 16), xlab=expression("Number of rows ("*10^3*")"), main="Column access", cex.axis=1)
dev.off()

incoming <- read.table("../timings_hdf5_row.txt", header=TRUE)
incoming$Ngenes <- incoming$Ngenes/1e3

pdf("HDF5_row_nrow.pdf")
plotter(incoming, "Ngenes", c("red", "black"), pch=c(18, 16), xlab=expression("Number of columns ("*10^3*")"), main="Row access", cex.axis=1, loc=NA)
dev.off()

# Layout type.

incoming <- read.table("../timings_hdf5_col_layout.txt", header=TRUE, sep="\t")
incoming <- incoming[-grep("uncompressed", incoming$Type),]
incoming$Type <- factor(incoming$Type, c("Contiguous", "Column chunks", "Row chunks", "Rechunked")) 

pdf("HDF5_col_layout.pdf")
plotter(incoming, "Ncells", c("grey50", "red", "tan4", "blue"), pch=c(16, 17, 15, 18), lty=c(1,1,2,2),
        xlab="Number of columns", main="Column access", cex.axis=1, upper=50000)
dev.off()

incoming <- read.table("../timings_hdf5_row_layout.txt", header=TRUE, sep="\t")
incoming <- incoming[-grep("uncompressed", incoming$Type),]
incoming$Type <- factor(incoming$Type, c("Contiguous", "Row chunks", "Column chunks", "Rechunked")) 
incoming$Time <- incoming$Time/1e3

pdf("HDF5_row_layout.pdf")
plotter(incoming, "Ngenes", c("grey50", "red", "tan4", "blue"), pch=c(16, 17, 15, 18), lty=c(1,1,2,2),
        xlab="Number of rows", main="Row access", cex.axis=1, upper=100, ylab="Time (s)", yaxt="n")
ticks <- c(10^(-2:2))
axis(2, at=ticks, ticks)
dev.off()


