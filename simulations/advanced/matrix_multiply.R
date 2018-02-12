library(beachtime)
library(Matrix)
library(HDF5Array)
overwrite <- TRUE

for (n in c(500, 1000, 2000)) {
    beach.dense.time <- default.dense.time <- 
        beach.sparse.time <- beach.sparse.time2 <- default.sparse.time <- 
        beach.hdf5.time <- beach.hdf5.time2 <- default.hdf5.time <- numeric(10)

    chunksize <- 100
    options(DelayedArray.block.size=n*chunksize*8) # Avoid just loading the entire matrix in.
    
    for (it in 1:10) {
        a <- matrix(runif(1e6), n, n)
        b <- matrix(runif(1e6), n, n)
        
        beach.dense.time[it] <- timeExprs(.Call(beachtime:::cxx_standard_matrix_multiplication, a, b))
        default.dense.time[it] <- timeExprs(a %*% b) 

        a3 <- rsparsematrix(1000, 1000, 0.01)
        b3 <- rsparsematrix(1000, 1000, 0.01)
        beach.sparse.time[it] <- timeExprs(.Call(beachtime:::cxx_standard_matrix_multiplication, a3, b3))
        beach.sparse.time2[it] <- timeExprs(.Call(beachtime:::cxx_indexed_matrix_multiplication, a3, b3))
        default.sparse.time[it] <- timeExprs(a3 %*% b3) 
    
        fpaths <- c("a.h5", "b.h5")
        a2 <- writeHDF5Array(a, fpaths[1], name="yay", chunk_dim=c(chunksize, chunksize), level=6)
        b2 <- writeHDF5Array(b, fpaths[2], name="yay", chunk_dim=c(chunksize, chunksize), level=6)
        beach.hdf5.time[it] <- timeExprs(.Call(beachtime:::cxx_standard_matrix_multiplication, a2, b), times=1) # Switching order is slower, as inner loop reads from file. 
        beach.hdf5.time2[it] <- timeExprs(.Call(beachtime:::cxx_standard_matrix_multiplication, a2, b2), times=1) 
        default.hdf5.time[it] <- timeExprs(a2 %*% b, times=1) # Multiplication of 2 DelayedArray objects is not supported.
        unlink(fpaths)
    }
    
    writeToFile(Type="dense (beachmat)", N=n, timings=beach.dense.time, file="timings_mat_mult.txt", overwrite=overwrite)
    overwrite <- FALSE 
    writeToFile(Type="dense (R)", N=n, timings=default.dense.time, file="timings_mat_mult.txt", overwrite=overwrite)
    writeToFile(Type="HDF5/dense (beachmat)", N=n, timings=beach.hdf5.time, file="timings_mat_mult.txt", overwrite=overwrite)
    writeToFile(Type="HDF5/HDF5 (beachmat)", N=n, timings=beach.hdf5.time2, file="timings_mat_mult.txt", overwrite=overwrite)
    writeToFile(Type="HDF5/dense (R)", N=n, timings=default.hdf5.time, file="timings_mat_mult.txt", overwrite=overwrite)
    writeToFile(Type="sparse (beachmat)", N=n, timings=beach.sparse.time, file="timings_mat_mult.txt", overwrite=overwrite)
    writeToFile(Type="sparse (beachmat II)", N=n, timings=beach.sparse.time2, file="timings_mat_mult.txt", overwrite=overwrite)
    writeToFile(Type="sparse (R)", N=n, timings=default.sparse.time, file="timings_mat_mult.txt", overwrite=overwrite)
}
