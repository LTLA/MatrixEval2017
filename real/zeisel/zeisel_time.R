# Timing on Zeisel's brain cell dataset.

infile <- "expression_mRNA_17-Aug-2014.txt"
counts <- read.delim(infile, stringsAsFactors=FALSE, header=FALSE, row.names=1, skip=11)[,-1] # First column after row names is some useless filler.

# Setting up the different representations:

dense.counts <- as.matrix(counts)
storage.mode(dense.counts) <- "double" # for proper comparison with sparse matrix.

library(Matrix)
sparse.counts <- as(Matrix(dense.counts), "dgCMatrix")

library(HDF5Array)
fpath <- c("bycol.h5", "byrow.h5", "rect.h5")
for (f in fpath) { if (file.exists(f)) { unlink(f) } }
hdf5.by.col <- writeHDF5Array(dense.counts, fpath[1], name="yay", chunk_dim=c(nrow(dense.counts), 1), level=6)
hdf5.by.row <- writeHDF5Array(dense.counts, fpath[2], name="yay", chunk_dim=c(1, ncol(dense.counts)), level=6)
hdf5.rect <- writeHDF5Array(dense.counts, fpath[3], name="yay", chunk_dim=c(200, 200), level=6)

write.table(data.frame(Sizes=c(dense=object.size(dense.counts)/1e3,
                               sparse=object.size(sparse.counts)/1e3,
                               HDF5.col=object.size(hdf5.by.col)/1e3,
                               HDF5.row=object.size(hdf5.by.row)/1e3,
                               HDF5.rect=object.size(hdf5.rect)/1e3,
                               HDF5.col.file=file.info(fpath[1])$size/1e3,
                               HDF5.row.file=file.info(fpath[2])$size/1e3,
                               HDF5.rect.file=file.info(fpath[3])$size/1e3)),                      ,
            col.names=NA, quote=FALSE, sep="\t", file="memory.txt")

# Running the column access timings.

library(beachtime)
dense.time <- sparse.time <- hdf5.time <- rect.time <- numeric(10)
for (it in 1:10) { 
    dense.time[it] <- timeExprs(BeachmatColSum(dense.counts), times=10)
    sparse.time[it] <- timeExprs(BeachmatColSum(sparse.counts), times=10)
    hdf5.time[it] <- timeExprs(BeachmatColSum(hdf5.by.col), times=1)
    rect.time[it] <- timeExprs(BeachmatColSum(hdf5.rect), times=1)
}

writeToFile(Type="dense", timings=dense.time, file="timings_col.txt", overwrite=TRUE)
writeToFile(Type="sparse", timings=sparse.time, file="timings_col.txt", overwrite=FALSE)
writeToFile(Type="HDF5 (column)", timings=hdf5.time, file="timings_col.txt", overwrite=FALSE)
writeToFile(Type="HDF5 (rectangle)", timings=rect.time, file="timings_col.txt", overwrite=FALSE)

# Running the row access timings.

dense.time <- sparse.time <- hdf5.time <- rect.time <- numeric(10)
for (it in 1:10) { 
    dense.time[it] <- timeExprs(BeachmatRowSum(dense.counts), times=10)
    sparse.time[it] <- timeExprs(BeachmatRowSum(sparse.counts), times=10)
    hdf5.time[it] <- timeExprs(BeachmatRowSum(hdf5.by.row), times=1)
    rect.time[it] <- timeExprs(BeachmatRowSum(hdf5.rect), times=1) 
}

writeToFile(Type="dense", timings=dense.time, file="timings_row.txt", overwrite=TRUE)
writeToFile(Type="sparse", timings=sparse.time, file="timings_row.txt", overwrite=FALSE)
writeToFile(Type="HDF5 (row)", timings=hdf5.time, file="timings_row.txt", overwrite=FALSE)
writeToFile(Type="HDF5 (rectangle)", timings=rect.time, file="timings_row.txt", overwrite=FALSE)

