# This makes bar plots of row/column access times.

pdf("zeisel_col.pdf")
par(mar=c(7.1, 5.1, 2.1, 2.1))
col.times <- read.table("../timings_col.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE)
barplot(setNames(col.times$Time, col.times$Type), ylab="Time (ms)", cex.axis=1, cex.lab=1.4, cex.names=1.4, 
        names.arg=sub(" ", "\n", col.times$Type), las=2,
        cex.main=1.4, main="Column access")
dev.off()

pdf("zeisel_row.pdf")
par(mar=c(7.1, 5.1, 2.1, 2.1))
row.times <- read.table("../timings_row.txt", header=TRUE, sep="\t")
barplot(setNames(row.times$Time, row.times$Type), ylab="Time (ms)", cex.axis=1, cex.lab=1.4, cex.names=1.4, 
        names.arg=sub(" ", "\n", row.times$Type), las=2,
        cex.main=1.4, main="Row access")
dev.off()
