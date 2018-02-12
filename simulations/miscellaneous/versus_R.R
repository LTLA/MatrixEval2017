library(beachtime)
library(Matrix)
library(HDF5Array)

block.sizes <- c(10, 20, 50, 100, 200, 500, 1000)

beach.time <- rep(list(numeric(10)), 3)
default.time <- rep(list(matrix(0, 10, length(block.sizes))), 3)
names(beach.time) <- names(default.time) <- c("dense", "sparse", "hdf5")

for (it in 1:10) { 
    for (mode in c("dense", "sparse", "hdf5")) { 
        if (mode=="dense") {
            X <- matrix(rnorm(1000000), nrow=100)
        } else if (mode=="sparse") {
            X <- rsparsematrix(nrow=100, ncol=10000, density=0.01)
        } else {
            fpath <- "a.h5"
            if (file.exists(fpath)) { unlink(fpath) }
            y <- matrix(rnorm(1000000), nrow=100)
            X <- writeHDF5Array(a, fpath, name="yay", chunk_dim=c(100, 100), level=6)
        }

        N <- ifelse(mode=="hdf5", 1, 10)
        beach.time[[mode]][it] <- timeExprs(BeachmatColSum(X), times=N)

        # Block processing with varying block sizes.
        # We use 'as.matrix' here to reflect the fact that we're using a 
        # single C(++) implementation of an arbitrary function, and 
        # coercing everything into a dense array to run it.
        for (b in seq_along(block.sizes)) { 
            block <- block.sizes[b]
            default.time[[mode]][it,b] <- timeExprs({
                nblocks <- ncol(X)/block
                for (x in seq_len(nblocks)) { 
                    base::colSums(as.matrix(X[,(x-1)*block + seq_len(block)]))
                }
            }, times=N)
        }
    }
}

overwrite <- TRUE
for(mode in names(beach.time)) {
    writeToFile(Type=mode, Block=NA, timings=beach.time[[mode]], file="timings_block.txt", overwrite=overwrite)
    overwrite <- FALSE
    for (b in seq_along(block.sizes)) {
        writeToFile(Type=mode, Block=block.sizes[b], timings=default.time[[mode]][,b], file="timings_block.txt", overwrite=overwrite)
    }
}

