library(beachtime)
tmp.dir <- "tmp-hdf5"
dir.create(tmp.dir)
library(HDF5Array)

############################
# This tests the optimality of the cache size calculation by repeating 
# runs after reducing the cache and recompiling. We do this for the simple
# case where the chunks are fully nested, as well as a case where the 
# chunks overshoot the edges of the matrix.

ngenes <- 10000    
ncells <- 2000
overwrite <- TRUE
fpaths <- file.path(tmp.dir, c("mult.h5", "nonmult.h5"))

for (redchunk in 0:1) {
    system(sprintf("g++ -DDECACHE=%i -std=c++11 -I%s -o HDF5ChunkTester cache_test.cpp %s -lz -ldl", 
           redchunk,
           system.file("include", package="Rhdf5lib"), 
           capture.output(Rhdf5lib::pkgconfig())))

    hdf5m.col.time <- hdf5m.row.time <- numeric(10)
    hdf5nm.col.time <- hdf5nm.row.time <- numeric(10)
    for (it in seq_len(10)) {         
        dense.counts <- matrix(rnorm(ngenes*ncells), ngenes, ncells)
        out.mult <- writeHDF5Array(dense.counts, fpaths[1], name="yay", chunkdim=c(100, 100), level=6)
        out.nmult <- writeHDF5Array(dense.counts, fpaths[2], name="yay", chunkdim=c(109, 97), level=6)

        hdf5m.col.time[it] <- timeExprs(system(sprintf("./HDF5ChunkTester %s yay 0", fpaths[1])), times=1)
        hdf5m.row.time[it] <- timeExprs(system(sprintf("./HDF5ChunkTester %s yay 1", fpaths[1])), times=1)
        hdf5nm.col.time[it] <- timeExprs(system(sprintf("./HDF5ChunkTester %s yay 0", fpaths[2])), times=1)
        hdf5nm.row.time[it] <- timeExprs(system(sprintf("./HDF5ChunkTester %s yay 1", fpaths[2])), times=1)
        unlink(fpaths)
    }

    writeToFile(Type=paste0("Nested:", redchunk), Ngenes=ngenes, Ncells=ncells,
                timings=hdf5m.row.time, file="timings_optim_cache_row.txt", overwrite=overwrite)
    writeToFile(Type=paste0("Nested:", redchunk), Ngenes=ngenes, Ncells=ncells,
                timings=hdf5m.col.time, file="timings_optim_cache_col.txt", overwrite=overwrite)
    overwrite <- FALSE 

    writeToFile(Type=paste0("Non-nested:", redchunk), Ngenes=ngenes, Ncells=ncells,
                timings=hdf5nm.row.time, file="timings_optim_cache_row.txt", overwrite=overwrite)
    writeToFile(Type=paste0("Non-nested:", redchunk), Ngenes=ngenes, Ncells=ncells,
                timings=hdf5nm.col.time, file="timings_optim_cache_col.txt", overwrite=overwrite)
}

############################
# This tests the optimality of the nslots calculation by repeating the  
# runs after replacing the nslots with the number of chunks per column.
# This shouldn't affect column access but row access should be totally broken.

for (redchunk in 0:1) {
    system(sprintf("g++ %s -std=c++11 -I%s -o HDF5ChunkTester tester.cpp %s -lz -ldl", 
           ifelse(redchunk==0, "", "-DBADSLOT"),
           system.file("include", package="Rhdf5lib"), 
           capture.output(Rhdf5lib::pkgconfig())))

    hdf5m.col.time <- hdf5m.row.time <- numeric(10)
    hdf5nm.col.time <- hdf5nm.row.time <- numeric(10)
    for (it in seq_len(10)) {         
        dense.counts <- matrix(rnorm(ngenes*ncells), ngenes, ncells)
        out.mult <- writeHDF5Array(dense.counts, fpaths[1], name="yay", chunkdim=c(100, 100), level=6)
        out.nmult <- writeHDF5Array(dense.counts, fpaths[2], name="yay", chunkdim=c(109, 97), level=6)

        hdf5m.col.time[it] <- timeExprs(system(sprintf("./HDF5ChunkTester %s yay 0", fpaths[1])), times=1)
        hdf5m.row.time[it] <- timeExprs(system(sprintf("./HDF5ChunkTester %s yay 1", fpaths[1])), times=1)
        hdf5nm.col.time[it] <- timeExprs(system(sprintf("./HDF5ChunkTester %s yay 0", fpaths[2])), times=1)
        hdf5nm.row.time[it] <- timeExprs(system(sprintf("./HDF5ChunkTester %s yay 1", fpaths[2])), times=1)
        unlink(fpaths)
    }

    writeToFile(Type=paste0("Nested:", redchunk), Ngenes=ngenes, Ncells=ncells,
                timings=hdf5m.row.time, file="timings_optim_nslot_row.txt", overwrite=overwrite)
    writeToFile(Type=paste0("Nested:", redchunk), Ngenes=ngenes, Ncells=ncells,
                timings=hdf5m.col.time, file="timings_optim_nslot_col.txt", overwrite=overwrite)
    overwrite <- FALSE 

    writeToFile(Type=paste0("Non-nested:", redchunk), Ngenes=ngenes, Ncells=ncells,
                timings=hdf5nm.row.time, file="timings_optim_nslot_row.txt", overwrite=overwrite)
    writeToFile(Type=paste0("Non-nested:", redchunk), Ngenes=ngenes, Ncells=ncells,
                timings=hdf5nm.col.time, file="timings_optim_nslot_col.txt", overwrite=overwrite)
}


