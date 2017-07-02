################################################################################
# Making a cell cycle picture.

assignments <- readRDS("../objects/cycle_output.rds")
X <- assignments$scores$G1
Y <- assignments$scores$G2M
X <- ceiling(X * 100)
Y <- ceiling(Y * 100)

output <- table(X, Y)
output <- log10(output+1)
my.cols <- rev(grey.colors(15))
upper <- ceiling(max(output))

png("cycle.png", width=8, height=7, units="in", res=300, pointsize=12)
layout(cbind(1,2), width=c(10, 2))
old.mar <- par()$mar
par(mar=c(5.1, 4.1, 4.1, 0.5))
image(output, col=my.cols, xlab="G1 score", ylab="G2M score", cex.axis=1.2, cex.lab=1.4, zlim=c(0, upper)) 
box()

segments(-1, 0.5, 0.5, 0.5, col="red", lty=2, lwd=2)
segments(0.5, -1, 0.5, 0.5, col="red", lty=2, lwd=2)
segments(0.5, 0.5, 2, 2, col="red", lty=2, lwd=2)

text(0, 0.6, col="red", sprintf("G2M phase:\n%i", sum(assignments$phase=="G2M")), cex=1.2, pos=4)
text(0, 0.4, col="red", sprintf("S phase:\n%i", sum(assignments$phase=="S")), cex=1.2, pos=4)
text(0.5, 0.4, col="red", sprintf("G1 phase:\n%i", sum(assignments$phase=="G1")), cex=1.2, pos=4)

# Adding the legend.
par(mar=c(5.1, 0, 4.1, 0))
plot.new()
plot.window(xlim=c(-0.1,2), ylim=c(0, 10))

start <- 4
inc <- 0.3
for (x in seq_along(my.cols)) { 
    rect(0, start+inc*(x-1), 0.5, start+inc*x, col=my.cols[x], border=my.cols[x])
}

top.height <- start+inc*length(my.cols)
text(0.5, start+inc/2, pos=4, 0)
text(0.5, top.height-inc/2, pos=4, upper)
text(0, top.height+inc, pos=4, expression(Log[10]*"(cells+1)"), cex=0.8, offset=0)
par(mar=old.mar)
dev.off()

################################################################################
# Making a plot of the size factors.

library(SummarizedExperiment)
se.out <- readRDS("../objects/qc_mat.rds")

ratio <- log(se.out$size_factors/se.out$Libsize)
nmads <- abs((ratio - median(ratio))/mad(ratio))

library(viridis)
my.cols <- rev(viridis(11))
coldex <- findInterval(nmads, 0:10/2)

options(bitmapType="cairo")
png("sizefacs.png", width=7, height=7, units="in", res=300, pointsize=12)
plot(se.out$Libsize/1e3, se.out$size_factors, log="xy", col=my.cols[coldex], 
     xlab=expression("Library size ("*10^3*")"), ylab="Size factor", 
     cex.axis=1.2, cex.lab=1.4, pch=16, cex=0.2)
dev.off()

