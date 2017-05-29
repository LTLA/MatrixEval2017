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
