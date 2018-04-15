library(beachtime)
library(Matrix)
library(HDF5Array)
library(DelayedMatrixStats)

infile <- "expression_mRNA_17-Aug-2014.txt"
counts <- read.delim(infile, stringsAsFactors=FALSE, header=FALSE, row.names=1, skip=11)[,-1] # First column after row names is some useless filler.

# Setting up the different representations:

dense.counts <- as.matrix(counts)
storage.mode(dense.counts) <- "double" # for proper comparison with sparse matrix.

sparse.counts <- as(Matrix(dense.counts), "dgCMatrix")

chunksize <- 200
hdf5.counts <- writeHDF5Array(dense.counts, chunkdim=c(chunksize, chunksize), level=6)
options(DelayedArray.block.size=nrow(counts)*chunksize*8) # Avoid just loading the entire matrix in, which would be misleading.

# Running across all modes.

for (mode in c("library_sizes", "detect_cells", "detect_genes")) { 
    if (mode=="library_sizes") {
        bFUN <- BeachmatColSumNoCopy 
        rFUNd <- function(M) { base::colSums(M) }
        rFUNs <- function(M) { Matrix::colSums(M) }
        rFUNh <- function(M) { DelayedMatrixStats::colSums2(M) }
    } else if (mode=="detect_cells") {
        bFUN <- numCellsExpressing
        rFUNd <- function(M) { base::rowSums(M > 0) }
        rFUNs <- function(M) { Matrix::rowSums(M > 0) }
        rFUNh <- function(M) { DelayedMatrixStats::rowSums2(M > 0) }
    } else {
        bFUN <- cellularDetectionRate
        rFUNd <- function(M) { base::colSums(M > 0) }
        rFUNs <- function(M) { Matrix::colSums(M > 0) }
        rFUNh <- function(M) { DelayedMatrixStats::colSums2(M > 0) }
    }

    beach.dense.time <- default.dense.time <- 
        beach.sparse.time <- default.sparse.time <- 
        beach.hdf5.time <- default.hdf5.time <- numeric(10)

    for (it in 1:10) {
        beach.dense.time[it] <- timeExprs(bFUN(dense.counts))
        default.dense.time[it] <- timeExprs(rFUNd(dense.counts))

        beach.sparse.time[it] <- timeExprs(bFUN(sparse.counts))
        default.sparse.time[it] <- timeExprs(rFUNs(sparse.counts)) 
    
        beach.hdf5.time[it] <- timeExprs(bFUN(hdf5.counts), times=1)
        default.hdf5.time[it] <- timeExprs(rFUNh(hdf5.counts), times=1)
    }
    
    overwrite <- TRUE
    fout <- paste0("timings_", mode, ".txt")
    writeToFile(Type="ordinary (beachmat)", timings=beach.dense.time, file=fout, overwrite=overwrite)
    overwrite <- FALSE 
    writeToFile(Type="ordinary (R)", timings=default.dense.time, file=fout, overwrite=overwrite)
    writeToFile(Type="HDF5 (beachmat)", timings=beach.hdf5.time, file=fout, overwrite=overwrite)
    writeToFile(Type="HDF5 (R)", timings=default.hdf5.time, file=fout, overwrite=overwrite)
    writeToFile(Type="sparse (beachmat)", timings=beach.sparse.time, file=fout, overwrite=overwrite)
    writeToFile(Type="sparse (R)", timings=default.sparse.time, file=fout, overwrite=overwrite)
}
