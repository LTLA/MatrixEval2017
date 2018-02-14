# This makes bar plots of row/column access times.

colors <- c(`ordinary`="darkblue",
            `sparse`="blue",
            `HDF5 (column)`="dodgerblue",
            `HDF5 (rectangle)`="lightblue")

pdf("zeisel_col.pdf")
par(mar=c(7.1, 5.1, 2.1, 2.1))
col.times <- read.table("../timings_col.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE)
barplot(setNames(col.times$Time/1000, col.times$Type), ylab="Time (s)", cex.axis=1, cex.lab=1.4, cex.names=1.4, 
        names.arg=sub(" ", "\n", col.times$Type), las=2, col=colors[col.times$Type],
        cex.main=1.4, main="Column access")
dev.off()

names(colors) <- sub("column", "row", names(colors))

pdf("zeisel_row.pdf")
par(mar=c(7.1, 5.1, 2.1, 2.1))
row.times <- read.table("../timings_row.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE)
barplot(setNames(row.times$Time/1000, row.times$Type), ylab="Time (s)", cex.axis=1, cex.lab=1.4, cex.names=1.4, 
        names.arg=sub(" ", "\n", row.times$Type), las=2, col=colors[row.times$Type],
        cex.main=1.4, main="Row access")
dev.off()

# This makes bar plots of timing comparisons to R.

colors <- c(`ordinary (beachmat)`="darkblue",
            `ordinary (R)`="darkblue",
            `sparse (beachmat)`="blue",
            `sparse (R)`="blue",
            `HDF5 (beachmat)`="dodgerblue",
            `HDF5 (R)`="dodgerblue")
pchs <- c(16, 4, 2)

pdf("zeisel_detection.pdf", width=10, height=6)
layout(cbind(1,2), width=c(5, 1.5))
par(mar=c(5.1, 4.1, 2.1, 0.1))
plot(1,1,type='n', xlim=c(1, 6.5), xlab="", xaxt="n", ylim=c(10, 2000), log="y", ylab="Time (ms)")
counter <- 0
for (mode in c("library_sizes", "detect_cells", "detect_genes")) { 
    incoming <- read.table(paste0("../timings_", mode, ".txt"), header=TRUE, sep="\t", stringsAsFactors=FALSE)
    vals <- setNames(incoming$Time, incoming$Type)
    vals <- vals[names(colors)]

    X <- seq_along(vals) + counter/5
    segments(X, vals, X, 1, lty=3)
    points(X, vals, col=colors, pch=pchs[counter+1], cex=1.5)
    counter <- counter + 1
}
axis(side=1, at=seq_along(colors)+1/5, labels=sub(" ", "\n", names(colors)), line=1, tick=FALSE)

par(mar=c(5.1, 0.1, 2.1, 0.1))
plot.new()
legend("left", col="black", pch=pchs, c("Library size per cell", "Number of cells per gene", "Number of genes per cell"), pt.cex=1.5)
dev.off()
