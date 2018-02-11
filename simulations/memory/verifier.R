# This checks that the memory usage for each matrix type is as expected.

library(Matrix)
library(HDF5Array)
overwrite <- TRUE

for (ngenes in c(1000, 2000, 5000, 10000)) {
    for (ncells in c(1000, 2000, 5000, 10000)) {

        dense.counts <- matrix(rnorm(ngenes*ncells), ngenes, ncells)
        write.table(data.frame(Ngenes=ngenes, Ncells=ncells, Size=as.numeric(object.size(dense.counts))/1e3),
                    file="memory_base.txt", append=!overwrite, col.names=overwrite, sep="\t", quote=FALSE, row.names=FALSE)

        hdf5.counts <- as(dense.counts, "HDF5Matrix")
        write.table(data.frame(Ngenes=ngenes, Ncells=ncells, Size=as.numeric(object.size(hdf5.counts))/1e3),
                    file="memory_hdf5.txt", append=!overwrite, col.names=overwrite, sep="\t", quote=FALSE, row.names=FALSE)
        write.table(data.frame(Ngenes=ngenes, Ncells=ncells, Size=file.info(hdf5.counts@seed@file)$size/1e3),
                    file="file_hdf5.txt", append=!overwrite, col.names=overwrite, sep="\t", quote=FALSE, row.names=FALSE)

        for (density in c(0.01, 0.05, 0.1, 0.2, 0.5)) {
            sparse.counts <- rsparsematrix(ngenes, ncells, density)
            write.table(data.frame(Ngenes=ngenes, Ncells=ncells, Density=density, Size=as.numeric(object.size(sparse.counts))/1e3),
                        file="memory_sparse.txt", append=!overwrite, col.names=overwrite, sep="\t", quote=FALSE, row.names=FALSE)
            overwrite <- FALSE 
        }

        gc()
    }
}

# Checking the effect of sparsity on the HDF5 file size.

library(beachtime)
tmp.dir <- "tmp-hdf5"
dir.create(tmp.dir, showWarnings=FALSE)

ncells <- 1000
ngenes <- 10000
overwrite <- TRUE
for (density in c(0.01, 0.05, 0.1, 0.2, 0.5, 1)) { 
    fpaths <- file.path(tmp.dir, c("contig.h5", "colcomp.h5", "colcomp0.h5"))
    colchunk.size <- colchunk0.size <- contig.size <- numeric(10)

    for (it in seq_len(10)) {         
        sparse.counts <- as.matrix(rsparsematrix(ngenes, ncells, density))
        out.contig <- writeHDF5Array(sparse.counts, fpaths[1], name="yay", level=0)
        out.bycol <- writeHDF5Array(sparse.counts, fpaths[2], name="yay", chunk_dim=c(ngenes, 1), level=6)
        out.bycol0 <- writeHDF5Array(sparse.counts, fpaths[3], name="yay", chunk_dim=c(ngenes, 1), level=0)

        contig.size[it] <- file.info(path(out.contig))$size/1e3
        colchunk.size[it] <- file.info(path(out.bycol))$size/1e3
        colchunk0.size[it] <- file.info(path(out.bycol0))$size/1e3

	    unlink(fpaths)
    }

    write.table(data.frame(Type="Contiguous", Density=density, Size=mean(contig.size), SE=se(contig.size)),
                file="file_hdf5_compress.txt", append=!overwrite, col.names=overwrite, sep="\t", quote=FALSE, row.names=FALSE)
    overwrite <- FALSE
    write.table(data.frame(Type="Column chunks", Density=density, Size=mean(colchunk.size), SE=se(colchunk.size)),
                file="file_hdf5_compress.txt", append=!overwrite, col.names=overwrite, sep="\t", quote=FALSE, row.names=FALSE)
    write.table(data.frame(Type="Column chunks, uncompressed", Density=density, Size=mean(colchunk0.size), SE=se(colchunk0.size)),
                file="file_hdf5_compress.txt", append=!overwrite, col.names=overwrite, sep="\t", quote=FALSE, row.names=FALSE)

}

unlink(tmp.dir, recursive=TRUE)
