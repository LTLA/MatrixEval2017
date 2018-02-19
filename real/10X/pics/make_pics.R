################################################################################
# Making a cell cycle picture.

assignments <- readRDS("../objects/cycle_output.rds")
X <- assignments$scores$G1
Y <- assignments$scores$G2M

pdf("cycle.pdf")
smoothScatter(X, Y, xlab="G1 score", ylab="G2M score", cex.axis=1.2, cex.lab=1.4)

segments(-1, 0.5, 0.5, 0.5, col="red", lty=2, lwd=2)
segments(0.5, -1, 0.5, 0.5, col="red", lty=2, lwd=2)
segments(0.5, 0.5, 2, 2, col="red", lty=2, lwd=2)

text(0, 0.6, col="red", sprintf("G2M phase:\n%i", sum(assignments$phase=="G2M")), cex=1.2, pos=4)
text(0, 0.4, col="red", sprintf("S phase:\n%i", sum(assignments$phase=="S")), cex=1.2, pos=4)
text(0.5, 0.4, col="red", sprintf("G1 phase:\n%i", sum(assignments$phase=="G1")), cex=1.2, pos=4)
dev.off()

################################################################################
# Making a plot of the size factors.

library(SingleCellExperiment)
sce <- readRDS("../objects/sce.rds")

ratio <- log(sizeFactors(sce)/sce$total_counts)
nmads <- abs((ratio - median(ratio))/mad(ratio))

library(viridis)
my.cols <- rev(viridis(11))
coldex <- findInterval(nmads, 0:10/2)

png("sizefacs_points.png", width=7, height=7, units="in", res=300, pointsize=12)
par(mar=c(0,0,0,0))
options(bitmapType="cairo")
plot(sce$total_counts/1e3, sizeFactors(sce), log="xy", col=my.cols[coldex], 
     xlab="", ylab="", pch=16, cex=0.5, axes=FALSE)
dev.off()

pdf("sizefacs.pdf")
plot(sce$total_counts/1e3, sizeFactors(sce), log="xy", col=my.cols[coldex], 
     xlab=expression("Library size ("*10^3*")"), ylab="Size factor", 
     cex.axis=1.2, cex.lab=1.4, pch=16, cex=0.2, type="n")

limits <- par("usr")
library(png)
img <- readPNG("sizefacs_points.png")
rasterImage(img, 10^limits[1], 10^limits[3], 10^limits[2], 10^limits[4])
box()
dev.off()

################################################################################
# Making a plot of the top HVGs.

hvg.out <- read.table("../objects/hvg_output_endog.txt", header=TRUE)
pois.out <- read.table("../objects/hvg_output_pois.txt", header=TRUE)
is.sig <- hvg.out$bio > 0

pdf("hvg.pdf")
par(mar=c(5.1, 5.1, 4.1, 2.1))
plot(hvg.out$mean, hvg.out$total, pch=16, col="grey", 
     xlab=expression("Mean log"[2]~"expression"),
     ylab=expression("Variance of log"[2]~"expression"), 
     cex.axis=1.2, cex.lab=1.4)
points(hvg.out$mean[is.sig], hvg.out$total[is.sig], col="orange", pch=16)

o <- order(hvg.out$mean)
lines(hvg.out$mean[o], hvg.out$tech[o], col="red", lwd=2)
o <- order(pois.out$mean)
lines(pois.out$mean[o], pois.out$tech[o], col="dodgerblue", lwd=2)

legend("topright", sprintf("%i HVGs out of %i genes", sum(is.sig & !is.na(is.sig)), length(is.sig)), bty="n", cex=1.1)

# Adding the top set of genes.
chosen <- hvg.out$Symbol[1:10]
xoffset <- c(0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2)
yoffset <- c(0., -0.1, 0.1, 0.1, 0.1, 0., 0.1, 0.1, 0.1, 0.)

for (i in seq_along(chosen)) {
    idex <- hvg.out$Symbol==chosen[i]
    xpos <- hvg.out$mean[idex]
    ypos <- hvg.out$total[idex]
    xpos2 <- xpos + xoffset[i]
    ypos2 <- ypos + yoffset[i]
    text(xpos2, ypos2, chosen[i], pos=4, offset=0.1, cex=1.1)
    segments(xpos, ypos, xpos2, ypos2)
}
dev.off()

################################################################################
# Making a PCA plot of the first two PCs.

#pc.data <- readRDS("../objects/sce.rds")
pc.mat <- reducedDim(sce)
pc.var <- attr(pc.mat, "percentVar")

pdf("pca.pdf")
smoothScatter(pc.mat[,1], pc.mat[,2], colramp=colorRampPalette(c("white", "black")),
     xlab=sprintf("PC1 (%.2f%%)", pc.var[1]*100),
     ylab=sprintf("PC2 (%.2f%%)", pc.var[2]*100),
     cex.axis=1.2, cex.lab=1.4)
dev.off()

