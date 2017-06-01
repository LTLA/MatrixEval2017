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
limit <- 500 
for (ncells in c(100, 200, 500, 1000, 2000, 5000, 10000)) {

    col.chunk.time <- col.rechunk.time <- col.xchunk.time <- col.contig.time <- def.time <- numeric(10)
    for (it in seq_len(10)) {         
        dense.counts <- matrix(rnorm(ngenes*ncells), ngenes, ncells)
        out.contig <- writeHDF5Array(dense.counts, file.path(tmp.dir, "contig.h5"), level=0)
        out.bycol <- writeHDF5Array(dense.counts, file.path(tmp.dir, "colcomp.h5"), chunk_dim=c(5000, 1), level=6)

        col.contig.time[it] <- timeExprs(BeachmatColSum(out.contig), times=1)
        col.chunk.time[it] <- timeExprs(BeachmatColSum(out.bycol), times=1)
        def.time[it] <- timeExprs(BeachmatColSum(dense.counts))

        if (ncells <= limit) { # Takes too long for a large number of cells, so skipping if that's the case.
            out.byrow <- writeHDF5Array(dense.counts, file.path(tmp.dir, "colrow.h5"), chunk_dim=c(1, min(ncells, 1000)), level=6)
            col.rechunk.time[it] <- timeExprs({
                out.col2 <- rechunk(out.byrow, file.path(tmp.dir, "colcomp2.h5"), byrow=FALSE)
                BeachmatColSum(out.col2)
            }, times=1)
            col.xchunk.time[it] <- timeExprs(BeachmatColSum(out.byrow), times=1)
        }
    }

    writeToFile(Type="HDF5 (contiguous)", Ngenes=ngenes, Ncells=ncells,
                timings=col.contig.time, file="timings_hdf5_col.txt", overwrite=overwrite)
    overwrite <- FALSE 
    writeToFile(Type="HDF5 (column-chunk)", Ngenes=ngenes, Ncells=ncells,
                timings=col.chunk.time, file="timings_hdf5_col.txt", overwrite=overwrite)
    writeToFile(Type="simple", Ngenes=ngenes, Ncells=ncells, 
                timings=def.time, file="timings_hdf5_col.txt", overwrite=overwrite)

    if (ncells <= limit) {
        writeToFile(Type="HDF5 (rechunked)", Ngenes=ngenes, Ncells=ncells,
                    timings=col.rechunk.time, file="timings_hdf5_col.txt", overwrite=overwrite)
        writeToFile(Type="HDF5 (row-chunk)", Ngenes=ngenes, Ncells=ncells,
                    timings=col.xchunk.time, file="timings_hdf5_col.txt", overwrite=overwrite)
    }
}

###########################
# Row access

ncells <- 1000
overwrite <- TRUE
limit <- 5000
for (ngenes in c(100, 200, 500, 10000, 20000, 50000, 100000)) {

    row.chunk.time <- row.xchunk.time <- row.contig.time <- row.rechunk.time <- def.time <- numeric(10)
    for (it in seq_len(10)) {         
        dense.counts <- matrix(rnorm(ngenes*ncells), ngenes, ncells)
        row.chunk.dim <- c(1, ncells)
        out.row <- writeHDF5Array(dense.counts, file.path(tmp.dir, "rowcomp.h5"), chunk_dim=row.chunk.dim, level=6)

        row.chunk.time[it] <- timeExprs(BeachmatRowSum(out.row), times=1)
        def.time[it] <- timeExprs(BeachmatRowSum(dense.counts))

        if (ngenes <= limit) { # Takes too long, just do this once.
            out.contig <- writeHDF5Array(dense.counts, file.path(tmp.dir, "contig.h5"), level=0)
            out.col <- writeHDF5Array(dense.counts, file.path(tmp.dir, "colcomp.h5"), chunk_dim=c(min(ngenes, 5000), 1), level=6)

            row.contig.time[it] <- timeExprs(BeachmatRowSum(out.contig), times=1)
            row.xchunk.time[it] <- timeExprs(BeachmatRowSum(out.bycol), times=1)
            row.rechunk.time[it] <- timeExprs({
                out.row2 <- rechunk(out.byrow, file.path(tmp.dir, "rowcomp2.h5"), byrow=TRUE)
                BeachmatRowSum(out.row2)
            }, times=1)
        }
    }

    writeToFile(Type="HDF5 (row-chunk)", Ngenes=ngenes, Ncells=ncells,
                timings=row.chunk.time, file="timings_hdf5_row.txt", overwrite=overwrite)
    overwrite <- FALSE 
    writeToFile(Type="Rcpp", Ngenes=ngenes, Ncells=ncells, 
                timings=def.time, file="timings_hdf5_row.txt", overwrite=overwrite)

    if (ngenes <= limit) {
        writeToFile(Type="HDF5 (rechunked)", Ngenes=ngenes, Ncells=ncells,
                    timings=row.rechunk.time, file="timings_hdf5_row.txt", overwrite=overwrite)
        writeToFile(Type="HDF5 (contiguous)", Ngenes=ngenes, Ncells=ncells,
                    timings=row.contig.time, file="timings_hdf5_row.txt", overwrite=overwrite)
        writeToFile(Type="HDF5 (column-chunk)", Ngenes=ngenes, Ncells=ncells,
                    timings=row.xchunk.time, file="timings_hdf5_row.txt", overwrite=overwrite)
    }
}

###########################
# Wrapping up

unlink(tmp.dir, recursive=TRUE)

