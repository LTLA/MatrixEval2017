library(beachtime)
library(beachmat)
tmp.dir <- "tmp-hdf5"
dir.create(tmp.dir)
library(HDF5Array)

# With respect to increasing numbers of columns.

ngenes <- 10000    
overwrite <- TRUE
for (ncells in c(1000, 2000, 5000, 10000)) {
    fpaths <- file.path(tmp.dir, c("colcomp.h5", "rowcomp.h5", "rowcomp2.h5", "colcomp2.h5"))

    hdf5.col.time <- hdf5.row.time <- numeric(10)
    for (it in seq_len(10)) {         
        dense.counts <- matrix(rnorm(ngenes*ncells), ngenes, ncells)
        out.col <- writeHDF5Array(dense.counts, fpaths[1], name="yay", chunk_dim=c(5000, 1), level=6)
        out.row <- writeHDF5Array(dense.counts, fpaths[2], name="yay", chunk_dim=c(1, min(ncells, 5000)), level=6)

        hdf5.col.time[it] <- timeExprs(rechunkByMargins(out.col, 5000, outfile=fpaths[3], byrow=TRUE), times=1)
        hdf5.row.time[it] <- timeExprs(rechunkByMargins(out.row, 5000, outfile=fpaths[4], byrow=FALSE), times=1)
        unlink(fpaths)
    }

    writeToFile(Type="Column to row", Ngenes=ngenes, Ncells=ncells,
                timings=hdf5.col.time, file="timings_rechunk_col.txt", overwrite=overwrite)
    overwrite <- FALSE 
    writeToFile(Type="Row to column", Ngenes=ngenes, Ncells=ncells,
                timings=hdf5.row.time, file="timings_rechunk_col.txt", overwrite=overwrite)

}

# With respect to increasing numbers of rows 

ncells <- 1000
overwrite <- TRUE
for (ngenes in c(10000, 20000, 50000, 100000)) {
    fpaths <- file.path(tmp.dir, c("colcomp.h5", "rowcomp.h5", "rowcomp2.h5", "colcomp2.h5"))

    hdf5.col.time <- hdf5.row.time <- numeric(10)
    for (it in seq_len(10)) {         
        dense.counts <- matrix(rnorm(ngenes*ncells), ngenes, ncells)
        out.col <- writeHDF5Array(dense.counts, fpaths[1], name="yay", chunk_dim=c(5000, 1), level=6)
        out.row <- writeHDF5Array(dense.counts, fpaths[2], name="yay", chunk_dim=c(1, min(ncells, 5000)), level=6)

        hdf5.col.time[it] <- timeExprs(rechunkByMargins(out.col, 5000, outfile=fpaths[3], byrow=TRUE), times=1)
        hdf5.row.time[it] <- timeExprs(rechunkByMargins(out.row, 5000, outfile=fpaths[4], byrow=FALSE), times=1)
        unlink(fpaths)
    }

    writeToFile(Type="Column to row", Ngenes=ngenes, Ncells=ncells,
                timings=hdf5.col.time, file="timings_rechunk_row.txt", overwrite=overwrite)
    overwrite <- FALSE 
    writeToFile(Type="Row to column", Ngenes=ngenes, Ncells=ncells,
                timings=hdf5.row.time, file="timings_rechunk_row.txt", overwrite=overwrite)

}

# Cleaning up

unlink(tmp.dir, recursive=TRUE)

########################################3
# Testing chunk cache settings: these two should give the same time,
# despite one of the chunk sizes not being a multiple of the other.
#
# library(beachmat); library(HDF5Array)
# ngenes <- 10000; ncells <- 1000
# dense.counts <- matrix(rnorm(ngenes*ncells), ngenes, ncells)
# basic <- writeHDF5Array(dense.counts, "tester.h5", name="yay", chunk_dim=c(100, 100), level=6)
# rewrite <- timeExprs(rechunkByMargins(basic, 100, "alpha.h5"), times=1)
# rewritex <- timeExprs(rechunkByMargins(basic, 101, "beta.h5"), times=1)

