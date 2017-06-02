# Simulating row and column access with a base matrix.

library(beachtime)
tmp.dir <- "tmp-hdf5"
dir.create(tmp.dir)
library(HDF5Array)

library(usingRhdf5lib)
rechunk <- function(incoming, outname, byrow=TRUE) {
    cxxfun <- ifelse(byrow, "pack_by_row", "pack_by_col")
    tolerance <- 1e8/8
    .Call(cxxfun, incoming@seed@file, incoming@seed@name, outname, tolerance, PACKAGE="usingRhdf5lib")
    HDF5Array(outname, incoming@seed@name)
}

###########################
# Column access

ngenes <- 10000    
overwrite <- TRUE
for (ncells in c(1000, 2000, 5000, 10000)) {
    fpath <- file.path(tmp.dir, "colcomp.h5")

    hdf5.time <- def.time <- numeric(10)
    for (it in seq_len(10)) {         
        dense.counts <- matrix(rnorm(ngenes*ncells), ngenes, ncells)
        out.col <- writeHDF5Array(dense.counts, fpath, name="yay", chunk_dim=c(5000, 1), level=6)
        hdf5.time[it] <- timeExprs(BeachmatColSum(out.col), times=1)
        def.time[it] <- timeExprs(BeachmatColSum(dense.counts))
        unlink(fpath)
    }

    writeToFile(Type="HDF5Array", Ngenes=ngenes, Ncells=ncells,
                timings=hdf5.time, file="timings_hdf5_col.txt", overwrite=overwrite)
    overwrite <- FALSE 
    writeToFile(Type="simple", Ngenes=ngenes, Ncells=ncells, 
                timings=def.time, file="timings_hdf5_col.txt", overwrite=overwrite)
}

###########################
# Row access

ncells <- 1000
overwrite <- TRUE
for (ngenes in c(10000, 20000, 50000, 100000)) {
    fpath <- file.path(tmp.dir, "rowcomp.h5")

    hdf5.time <- def.time <- numeric(10)
    for (it in seq_len(10)) {         
        dense.counts <- matrix(rnorm(ngenes*ncells), ngenes, ncells)
        out.row <- writeHDF5Array(dense.counts, fpath, name="yay", chunk_dim=c(1, ncells), level=6)
        hdf5.time[it] <- timeExprs(BeachmatRowSum(out.row), times=1)
        def.time[it] <- timeExprs(BeachmatRowSum(dense.counts))
        unlink(fpath)
    }

    writeToFile(Type="HDF5Array", Ngenes=ngenes, Ncells=ncells,
                timings=hdf5.time, file="timings_hdf5_row.txt", overwrite=overwrite)
    overwrite <- FALSE 
    writeToFile(Type="simple", Ngenes=ngenes, Ncells=ncells, 
                timings=def.time, file="timings_hdf5_row.txt", overwrite=overwrite)
}

###########################
# Column access, alternative layouts.

ngenes <- 1000
overwrite <- TRUE
for (ncells in c(100, 200, 500)) { # a _lot_ shorter, as it takes too long using the full set of rows/columns.
    fpaths <- file.path(tmp.dir, c("contig.h5", "colcomp.h5", "rowcomp.h5", "rechunk.h5"))
    colchunk.time <- rowchunk.time <- contig.time <- rechunk.time <- numeric(10)

    for (it in seq_len(10)) {         
        dense.counts <- matrix(rnorm(ngenes*ncells), ngenes, ncells)
        out.contig <- writeHDF5Array(dense.counts, fpaths[1], name="yay", level=0)
        out.bycol <- writeHDF5Array(dense.counts, fpaths[2], name="yay", chunk_dim=c(ngenes, 1), level=6)
        out.byrow <- writeHDF5Array(dense.counts, fpaths[3], name="yay", chunk_dim=c(1, ncells), level=6)

        contig.time[it] <- timeExprs(BeachmatColSum(out.contig))
        colchunk.time[it] <- timeExprs(BeachmatColSum(out.bycol))
        rowchunk.time[it] <- timeExprs(BeachmatColSum(out.byrow), times=1)
        rechunk.time[it] <- timeExprs({
            out.col2 <- rechunk(out.byrow, fpaths[4], byrow=FALSE)
            BeachmatColSum(out.col2)
        }, times=1)

        unlink(fpaths)
    }

    writeToFile(Type="HDF5 (contiguous)", Ngenes=ngenes, Ncells=ncells,
                timings=contig.time, file="timings_hdf5_col_layout.txt", overwrite=overwrite)
    overwrite <- FALSE 
    writeToFile(Type="HDF5 (column-chunk)", Ngenes=ngenes, Ncells=ncells,
                timings=colchunk.time, file="timings_hdf5_col_layout.txt", overwrite=overwrite)
    writeToFile(Type="HDF5 (row-chunk)", Ngenes=ngenes, Ncells=ncells,
                timings=rowchunk.time, file="timings_hdf5_col_layout.txt", overwrite=overwrite)
    writeToFile(Type="HDF5 (rechunked)", Ngenes=ngenes, Ncells=ncells,
                timings=rechunk.time, file="timings_hdf5_col_layout.txt", overwrite=overwrite)
}

###########################
# Row access, alternative layouts.

ncells <- 100
overwrite <- TRUE
for (ngenes in c(1000, 2000, 5000)) { 
    fpaths <- file.path(tmp.dir, c("contig.h5", "colcomp.h5", "rowcomp.h5", "rechunk.h5"))
    colchunk.time <- rowchunk.time <- contig.time <- rechunk.time <- numeric(10)

    for (it in seq_len(10)) {
        dense.counts <- matrix(rnorm(ngenes*ncells), ngenes, ncells)
        out.contig <- writeHDF5Array(dense.counts, fpaths[1], name="yay", level=0)
        out.bycol <- writeHDF5Array(dense.counts, fpaths[2], name="yay", chunk_dim=c(ngenes, 1), level=6)
        out.byrow <- writeHDF5Array(dense.counts, fpaths[3], name="yay", chunk_dim=c(1, ncells), level=6) 
    
        contig.time[it] <- timeExprs(BeachmatRowSum(out.contig), times=1)
        colchunk.time[it] <- timeExprs(BeachmatRowSum(out.bycol), times=1)
        rowchunk.time[it] <- timeExprs(BeachmatRowSum(out.byrow), times=1)
        rechunk.time[it] <- timeExprs({
            out.row2 <- rechunk(out.bycol, fpaths[4], byrow=TRUE)
            BeachmatRowSum(out.row2)
        }, times=1)

        unlink(fpaths)
    }

    writeToFile(Type="HDF5 (contiguous)", Ngenes=ngenes, Ncells=ncells,
                timings=contig.time, file="timings_hdf5_row_layout.txt", overwrite=overwrite)
    overwrite <- FALSE 
    writeToFile(Type="HDF5 (column-chunk)", Ngenes=ngenes, Ncells=ncells,
                timings=colchunk.time, file="timings_hdf5_row_layout.txt", overwrite=overwrite)
    writeToFile(Type="HDF5 (row-chunk)", Ngenes=ngenes, Ncells=ncells,
                timings=rowchunk.time, file="timings_hdf5_row_layout.txt", overwrite=overwrite)
    writeToFile(Type="HDF5 (rechunked)", Ngenes=ngenes, Ncells=ncells,
                timings=rechunk.time, file="timings_hdf5_row_layout.txt", overwrite=overwrite)
}

###########################
# Wrapping up

unlink(tmp.dir, recursive=TRUE)

