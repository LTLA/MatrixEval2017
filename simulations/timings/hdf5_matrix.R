# Simulating row and column access with a base matrix.

library(beachtime)
tmp.dir <- "tmp-hdf5"
dir.create(tmp.dir)
library(HDF5Array)

makeContiguous <- function(data, file, name) {
    # Rewriting this function, as h5createDataset2 automatically sets the chunk dimensions to 'dim'.
    # I think the chunk dimensions are ignored by h5createDataset if level=0, but this just makes sure
    # by explicitly specifying chunk dimensions of NULL in h5createDataset().
    dim <- dim(data)
    type <- storage.mode(data)
    chunk_dim <- NULL
    level <- 0
    h5createFile(file)
    h5createDataset(file, name, dim, storage.mode=type, size=NULL, chunk = chunk_dim, level=level)
    appendDatasetCreationToHDF5DumpLog(file, name, dim, type, chunk_dim, level)
    h5write(data, file, name)
    HDF5Array(file, name)
}

###########################
# Column access

ngenes <- 10000    
overwrite <- TRUE
for (ncells in c(1000, 2000, 5000, 10000)) {
    fpaths <- file.path(tmp.dir, c("colcomp.h5", "sqcomp.h5"))

    hdf5.col.time <- hdf5.sq.time <- def.time <- numeric(10)
    for (it in seq_len(10)) {         
        dense.counts <- matrix(rnorm(ngenes*ncells), ngenes, ncells)
        out.col <- writeHDF5Array(dense.counts, fpaths[1], name="yay", chunk_dim=c(ngenes, 1), level=6)
        out.sq <- writeHDF5Array(dense.counts, fpaths[2], name="yay", chunk_dim=c(100, 100), level=6)

        hdf5.col.time[it] <- timeExprs(BeachmatColSum(out.col), times=1)
        hdf5.sq.time[it] <- timeExprs(BeachmatColSum(out.sq), times=1)
        def.time[it] <- timeExprs(BeachmatColSum(dense.counts))
        unlink(fpaths)
    }

    writeToFile(Type="HDF5 (column)", Ngenes=ngenes, Ncells=ncells,
                timings=hdf5.col.time, file="timings_hdf5_col.txt", overwrite=overwrite)
    overwrite <- FALSE 
    writeToFile(Type="HDF5 (rectangle)", Ngenes=ngenes, Ncells=ncells,
                timings=hdf5.sq.time, file="timings_hdf5_col.txt", overwrite=overwrite)
    writeToFile(Type="simple", Ngenes=ngenes, Ncells=ncells, 
                timings=def.time, file="timings_hdf5_col.txt", overwrite=overwrite)
}

###########################
# Row access

ncells <- 1000
overwrite <- TRUE
for (ngenes in c(10000, 20000, 50000, 100000)) {
    fpath <- file.path(tmp.dir, c("rowcomp.h5", "rect.h5"))

    hdf5.row.time <- hdf5.rect.time <- def.time <- numeric(10)
    for (it in seq_len(10)) {         
        dense.counts <- matrix(rnorm(ngenes*ncells), ngenes, ncells)
        out.row <- writeHDF5Array(dense.counts, fpaths[1], name="yay", chunk_dim=c(1, ncells), level=6)
        out.rect <- writeHDF5Array(dense.counts, fpaths[2], name="yay", chunk_dim=c(100, 100), level=6)

        hdf5.row.time[it] <- timeExprs(BeachmatRowSum(out.row), times=1)
        hdf5.rect.time[it] <- timeExprs(BeachmatRowSum(out.rect), times=1)
        def.time[it] <- timeExprs(BeachmatRowSum(dense.counts))
        unlink(fpaths)
    }

    writeToFile(Type="HDF5 (row)", Ngenes=ngenes, Ncells=ncells,
                timings=hdf5.row.time, file="timings_hdf5_row.txt", overwrite=overwrite)
    overwrite <- FALSE 
    writeToFile(Type="HDF5 (rectangle)", Ngenes=ngenes, Ncells=ncells,
                timings=hdf5.rect.time, file="timings_hdf5_row.txt", overwrite=overwrite)
    writeToFile(Type="simple", Ngenes=ngenes, Ncells=ncells, 
                timings=def.time, file="timings_hdf5_row.txt", overwrite=overwrite)
}

###########################
# Column access, alternative layouts.

ngenes <- 1000
overwrite <- TRUE
for (ncells in c(100, 200, 500)) { # a _lot_ shorter, as it takes too long using the full set of rows/columns.
    fpaths <- file.path(tmp.dir, c("contig.h5", "colcomp.h5", "colcomp0.h5", "rowcomp.h5", "rect.h5"))
    colchunk.time <- colchunk0.time <- rowchunk.time <- contig.time <- rect.time <- numeric(10)

    for (it in seq_len(10)) {         
        dense.counts <- matrix(rnorm(ngenes*ncells), ngenes, ncells)
        out.contig <- makeContiguous(dense.counts, fpaths[1], name="yay")
        out.bycol <- writeHDF5Array(dense.counts, fpaths[2], name="yay", chunk_dim=c(ngenes, 1), level=6)
        out.bycol0 <- writeHDF5Array(dense.counts, fpaths[3], name="yay", chunk_dim=c(ngenes, 1), level=1) # Reduced compression (level=0 doesn't actually chunk).
        out.byrow <- writeHDF5Array(dense.counts, fpaths[4], name="yay", chunk_dim=c(1, ncells), level=6)
        out.rect <- writeHDF5Array(dense.counts, fpaths[5], name="yay", chunk_dim=c(40, 40), level=6)

        contig.time[it] <- timeExprs(BeachmatColSum(out.contig))
        colchunk.time[it] <- timeExprs(BeachmatColSum(out.bycol))
        colchunk0.time[it] <- timeExprs(BeachmatColSum(out.bycol0))
        rowchunk.time[it] <- timeExprs(BeachmatColSum(out.byrow), times=1)
        rect.time[it] <- timeExprs(BeachmatColSum(out.rect), times=1)

        unlink(fpaths)
    }

    writeToFile(Type="Contiguous", Ngenes=ngenes, Ncells=ncells,
                timings=contig.time, file="timings_hdf5_col_layout.txt", overwrite=overwrite)
    overwrite <- FALSE 
    writeToFile(Type="Column chunks", Ngenes=ngenes, Ncells=ncells,
                timings=colchunk.time, file="timings_hdf5_col_layout.txt", overwrite=overwrite)
    writeToFile(Type="Column chunks, uncompressed", Ngenes=ngenes, Ncells=ncells,
                timings=colchunk0.time, file="timings_hdf5_col_layout.txt", overwrite=overwrite)
    writeToFile(Type="Row chunks", Ngenes=ngenes, Ncells=ncells,
                timings=rowchunk.time, file="timings_hdf5_col_layout.txt", overwrite=overwrite)
    writeToFile(Type="Rectangular chunks", Ngenes=ngenes, Ncells=ncells,
                timings=rect.time, file="timings_hdf5_col_layout.txt", overwrite=overwrite)
}

###########################
# Row access, alternative layouts.

ncells <- 100
overwrite <- TRUE
for (ngenes in c(1000, 2000, 5000)) { 
    fpaths <- file.path(tmp.dir, c("contig.h5", "colcomp.h5", "rowcomp.h5", "rowcomp0.h5", "rect.h5"))
    colchunk.time <- rowchunk.time <- rowchunk0.time <- contig.time <- rect.time <- numeric(10)

    for (it in seq_len(10)) {
        dense.counts <- matrix(rnorm(ngenes*ncells), ngenes, ncells)
        out.contig <- makeContiguous(dense.counts, fpaths[1], name="yay")
        out.bycol <- writeHDF5Array(dense.counts, fpaths[2], name="yay", chunk_dim=c(ngenes, 1), level=6)
        out.byrow <- writeHDF5Array(dense.counts, fpaths[3], name="yay", chunk_dim=c(1, ncells), level=6) 
        out.byrow0 <- writeHDF5Array(dense.counts, fpaths[4], name="yay", chunk_dim=c(1, ncells), level=1) # Reduced compression.
        out.rect <- writeHDF5Array(dense.counts, fpaths[5], name="yay", chunk_dim=c(40, 40), level=6)
 
        contig.time[it] <- timeExprs(BeachmatRowSum(out.contig), times=1)
        colchunk.time[it] <- timeExprs(BeachmatRowSum(out.bycol), times=1)
        rowchunk.time[it] <- timeExprs(BeachmatRowSum(out.byrow), times=1)
        rowchunk0.time[it] <- timeExprs(BeachmatRowSum(out.byrow0), times=1)
        rect.time[it] <- timeExprs(BeachmatRowSum(out.rect), times=1)

        unlink(fpaths)
    }

    writeToFile(Type="Contiguous", Ngenes=ngenes, Ncells=ncells,
                timings=contig.time, file="timings_hdf5_row_layout.txt", overwrite=overwrite)
    overwrite <- FALSE 
    writeToFile(Type="Column chunks", Ngenes=ngenes, Ncells=ncells,
                timings=colchunk.time, file="timings_hdf5_row_layout.txt", overwrite=overwrite)
    writeToFile(Type="Row chunks", Ngenes=ngenes, Ncells=ncells,
                timings=rowchunk.time, file="timings_hdf5_row_layout.txt", overwrite=overwrite)
    writeToFile(Type="Row chunks, uncompressed", Ngenes=ngenes, Ncells=ncells,
                timings=rowchunk0.time, file="timings_hdf5_row_layout.txt", overwrite=overwrite)
    writeToFile(Type="Rectangular chunks", Ngenes=ngenes, Ncells=ncells,
                timings=rect.time, file="timings_hdf5_row_layout.txt", overwrite=overwrite)
}

###########################
# Random column access, alternative layouts.

ngenes <- 1000
overwrite <- TRUE
for (ncells in c(100, 200, 500)) { # a _lot_ shorter, as it takes too long using the full set of rows/columns.
    fpaths <- file.path(tmp.dir, c("contig.h5", "colcomp.h5", "rowcomp.h5", "rect.h5"))
    colchunk.time <- rowchunk.time <- contig.time <- rect.time <- numeric(10)

    for (it in seq_len(10)) {         
        dense.counts <- matrix(rnorm(ngenes*ncells), ngenes, ncells)
        out.contig <- makeContiguous(dense.counts, fpaths[1], name="yay")
        out.bycol <- writeHDF5Array(dense.counts, fpaths[2], name="yay", chunk_dim=c(ngenes, 1), level=6)
        out.byrow <- writeHDF5Array(dense.counts, fpaths[3], name="yay", chunk_dim=c(1, ncells), level=6)
        out.rect <- writeHDF5Array(dense.counts, fpaths[4], name="yay", chunk_dim=c(40, 40), level=6)

        o <- sample(ncells)
        contig.time[it] <- timeExprs(BeachmatColSumRandom(out.contig, o))
        colchunk.time[it] <- timeExprs(BeachmatColSumRandom(out.bycol, o))
        rowchunk.time[it] <- timeExprs(BeachmatColSumRandom(out.byrow, o), times=1)
        rect.time[it] <- timeExprs(BeachmatColSumRandom(out.rect, o), times=1)

        unlink(fpaths)
    }

    writeToFile(Type="Contiguous", Ngenes=ngenes, Ncells=ncells,
                timings=contig.time, file="timings_hdf5_col_layout_random.txt", overwrite=overwrite)
    overwrite <- FALSE 
    writeToFile(Type="Column chunks", Ngenes=ngenes, Ncells=ncells,
                timings=colchunk.time, file="timings_hdf5_col_layout_random.txt", overwrite=overwrite)
    writeToFile(Type="Row chunks", Ngenes=ngenes, Ncells=ncells,
                timings=rowchunk.time, file="timings_hdf5_col_layout_random.txt", overwrite=overwrite)
    writeToFile(Type="Rectangular chunks", Ngenes=ngenes, Ncells=ncells,
                timings=rect.time, file="timings_hdf5_col_layout_random.txt", overwrite=overwrite)
}

###########################
# Random row access, alternative layouts.

ncells <- 100
overwrite <- TRUE
for (ngenes in c(1000, 2000, 5000)) { 
    fpaths <- file.path(tmp.dir, c("contig.h5", "colcomp.h5", "rowcomp.h5", "rect.h5"))
    colchunk.time <- rowchunk.time <- contig.time <- rect.time <- numeric(10)

    for (it in seq_len(10)) {
        dense.counts <- matrix(rnorm(ngenes*ncells), ngenes, ncells)
        out.contig <- makeContiguous(dense.counts, fpaths[1], name="yay")
        out.bycol <- writeHDF5Array(dense.counts, fpaths[2], name="yay", chunk_dim=c(ngenes, 1), level=6)
        out.byrow <- writeHDF5Array(dense.counts, fpaths[3], name="yay", chunk_dim=c(1, ncells), level=6) 
        out.rect <- writeHDF5Array(dense.counts, fpaths[4], name="yay", chunk_dim=c(40, 40), level=6)
 
        o <- order(ngenes)
        contig.time[it] <- timeExprs(BeachmatRowSumRandom(out.contig, o), times=1)
        colchunk.time[it] <- timeExprs(BeachmatRowSumRandom(out.bycol, o), times=1)
        rowchunk.time[it] <- timeExprs(BeachmatRowSumRandom(out.byrow, o), times=1)
        rect.time[it] <- timeExprs(BeachmatRowSumRandom(out.rect, o), times=1)

        unlink(fpaths)
    }

    writeToFile(Type="Contiguous", Ngenes=ngenes, Ncells=ncells,
                timings=contig.time, file="timings_hdf5_row_layout_random.txt", overwrite=overwrite)
    overwrite <- FALSE 
    writeToFile(Type="Column chunks", Ngenes=ngenes, Ncells=ncells,
                timings=colchunk.time, file="timings_hdf5_row_layout_random.txt", overwrite=overwrite)
    writeToFile(Type="Row chunks", Ngenes=ngenes, Ncells=ncells,
                timings=rowchunk.time, file="timings_hdf5_row_layout_random.txt", overwrite=overwrite)
    writeToFile(Type="Rectangular chunks", Ngenes=ngenes, Ncells=ncells,
                timings=rect.time, file="timings_hdf5_row_layout_random.txt", overwrite=overwrite)
}

###########################
# Testing the effects of sparsity on speed of column access.

library(Matrix)
ncells <- 1000
ngenes <- 10000
overwrite <- TRUE
for (density in c(0.01, 0.05, 0.1, 0.2, 0.5, 1)) { 
    fpaths <- file.path(tmp.dir, c("contig.h5", "colcomp.h5", "colcomp0.h5"))
    colchunk.time <- colchunk0.time <- contig.time <- numeric(10)

    for (it in seq_len(10)) {         
        sparse.counts <- as.matrix(rsparsematrix(ngenes, ncells, density))
        out.contig <- makeContiguous(sparse.counts, fpaths[1], name="yay")
        out.bycol <- writeHDF5Array(sparse.counts, fpaths[2], name="yay", chunk_dim=c(ngenes, 1), level=6)
        out.bycol0 <- writeHDF5Array(sparse.counts, fpaths[3], name="yay", chunk_dim=c(ngenes, 1), level=1) # Reduced compression.

        contig.time[it] <- timeExprs(BeachmatColSum(out.contig))
        colchunk.time[it] <- timeExprs(BeachmatColSum(out.bycol))
        colchunk0.time[it] <- timeExprs(BeachmatColSum(out.bycol0))

    	unlink(fpaths)
    }
	
    writeToFile(Type="Contiguous", Ngenes=ngenes, Ncells=ncells, Density=density,
                timings=contig.time, file="timings_hdf5_col_compress.txt", overwrite=overwrite)
    overwrite <- FALSE 
    writeToFile(Type="Column chunks", Ngenes=ngenes, Ncells=ncells, Density=density,
                timings=colchunk.time, file="timings_hdf5_col_compress.txt", overwrite=overwrite)
    writeToFile(Type="Column chunks, uncompressed", Ngenes=ngenes, Ncells=ncells, Density=density,
                timings=colchunk0.time, file="timings_hdf5_col_compress.txt", overwrite=overwrite)
}
 
###########################
# Wrapping up

unlink(tmp.dir, recursive=TRUE)
sessionInfo()
