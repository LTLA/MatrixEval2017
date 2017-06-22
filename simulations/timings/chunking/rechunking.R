library(beachtime)
library(scater)
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

        hdf5.col.time[it] <- timeExprs(scater:::rechunk(out.col, 5000, outfile=fpaths[3], byrow=TRUE), times=1)
        hdf5.row.time[it] <- timeExprs(scater:::rechunk(out.row, 5000, outfile=fpaths[4], byrow=FALSE), times=1)
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

        hdf5.col.time[it] <- timeExprs(scater:::rechunk(out.col, 5000, outfile=fpaths[3], byrow=TRUE), times=1)
        hdf5.row.time[it] <- timeExprs(scater:::rechunk(out.row, 5000, outfile=fpaths[4], byrow=FALSE), times=1)
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
