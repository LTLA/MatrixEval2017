# This makes bar plots of row/column access times.

pdf("zeisel_col.pdf")
col.times <- read.table("../timings_col.txt", header=TRUE)
barplot(setNames(col.times$Time, col.times$Type), ylab="Time (ms)", cex.axis=1, cex.lab=1.4, cex.names=1.4, 
        cex.main=1.4, main="Column access", log="y", ylim=range(col.times$Time)*c(0.5, 1.5))
dev.off()

pdf("zeisel_row.pdf")
row.times <- read.table("../timings_row.txt", header=TRUE)
barplot(setNames(row.times$Time, row.times$Type), ylab="Time (ms)", cex.axis=1, cex.lab=1.4, cex.names=1.4, 
        cex.main=1.4, main="Row access", log="y", ylim=range(row.times$Time)*c(0.5, 1.5))
dev.off()
