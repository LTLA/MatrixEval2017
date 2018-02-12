library(beachtime)
library(Matrix)
library(HDF5Array)
library(DelayedMatrixStats)

chunksize <- 100
Ngenes <- 10000
Ncells <- 1000
options(DelayedArray.block.size=Ngenes*chunksize*8) # Avoid just loading the entire matrix in, which would be misleading.

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
        a <- matrix(runif(1e6), Ngenes, Ncells)
        beach.dense.time[it] <- timeExprs(bFUN(a))
        default.dense.time[it] <- timeExprs(rFUNd(a))

        a3 <- rsparsematrix(Ngenes, Ncells, 0.01)
        beach.sparse.time[it] <- timeExprs(bFUN(a3))
        default.sparse.time[it] <- timeExprs(rFUNs(a3)) 
    
        fpaths <- "det_tmp.h5"
        if (file.exists(fpaths)) { unlink(fpaths) }
        a2 <- writeHDF5Array(a, fpaths[1], name="yay", chunk_dim=c(chunksize, chunksize), level=6)
        beach.hdf5.time[it] <- timeExprs(bFUN(a2), times=1)
        default.hdf5.time[it] <- timeExprs(rFUNh(a2), times=1)
        unlink(fpaths) 
    }
    
    overwrite <- TRUE
    fout <- paste0("timings_", mode, ".txt")
    writeToFile(Type="dense (beachmat)", timings=beach.dense.time, file=fout, overwrite=overwrite)
    overwrite <- FALSE 
    writeToFile(Type="dense (R)", timings=default.dense.time, file=fout, overwrite=overwrite)
    writeToFile(Type="HDF5 (beachmat)", timings=beach.hdf5.time, file=fout, overwrite=overwrite)
    writeToFile(Type="HDF5 (R)", timings=default.hdf5.time, file=fout, overwrite=overwrite)
    writeToFile(Type="sparse (beachmat)", timings=beach.sparse.time, file=fout, overwrite=overwrite)
    writeToFile(Type="sparse (R)", timings=default.sparse.time, file=fout, overwrite=overwrite)
}
