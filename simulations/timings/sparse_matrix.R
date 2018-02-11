# Simulating row and column access with a sparse matrix.

library(beachtime)
library(Matrix)

###########################
# Consecutive column access

ngenes <- 10000    
overwrite <- TRUE
for (ncells in c(1000, 2000, 5000, 10000)) {
    for (density in c(0.01, 0.05, 0.1, 0.2)) { 
        
        col.time <- colnocopy.time <- arma.time <- eigen.time <- def.time <- numeric(10)
        for (it in seq_len(10)) { 
            sparse.counts <- rsparsematrix(ngenes, ncells, density)
            dense.counts <- as.matrix(sparse.counts)
            col.time[it] <- timeExprs(BeachmatColSum(sparse.counts))
            colnocopy.time[it] <- timeExprs(BeachmatColSumNoCopy(sparse.counts))
            arma.time[it] <- timeExprs(ArmaColSum(sparse.counts))
            eigen.time[it] <- timeExprs(ArmaColSum(sparse.counts))
            def.time[it] <- timeExprs(BeachmatColSum(dense.counts))
        }
    
        writeToFile(Type="sparse", Ngenes=ngenes, Ncells=ncells, Density=density, 
                    timings=col.time, file="timings_sparse_col.txt", overwrite=overwrite)
        overwrite <- FALSE 
        writeToFile(Type="sparse (no copy)", Ngenes=ngenes, Ncells=ncells, Density=density, 
                    timings=colnocopy.time, file="timings_sparse_col.txt", overwrite=overwrite)
        writeToFile(Type="RcppArmadillo", Ngenes=ngenes, Ncells=ncells, Density=density, 
                    timings=arma.time, file="timings_sparse_col.txt", overwrite=overwrite)
        writeToFile(Type="RcppEigen", Ngenes=ngenes, Ncells=ncells, Density=density, 
                    timings=eigen.time, file="timings_sparse_col.txt", overwrite=overwrite)
        writeToFile(Type="simple", Ngenes=ngenes, Ncells=ncells, Density=density, 
                    timings=def.time, file="timings_sparse_col.txt", overwrite=overwrite)
    }
}

###########################
# Consecutive row access

ncells <- 1000
overwrite <- TRUE
for (ngenes in c(10000, 20000, 50000, 100000)) {
    for (density in c(0.01, 0.05, 0.1, 0.2)) { 
        skip.arma <- !(density==0.01 || ngenes<=50000)

        naive.time <- row.time <- arma.time <- eigen.time <- def.time <- numeric(10)
        for (it in seq_len(10)) { 
            sparse.counts <- rsparsematrix(ngenes, ncells, density)
            dense.counts <- as.matrix(sparse.counts)
            row.time[it] <- timeExprs(BeachmatRowSum(sparse.counts))
            def.time[it] <- timeExprs(BeachmatRowSum(dense.counts))
            naive.time[it] <- timeExprs(NaiveSparseRowSum(sparse.counts), times=1)

            # This takes a lot of time, so will only run once. 
            # We also skip the more complex timings. 
            if (!skip.arma) {
                eigen.time[it] <- timeExprs(EigenRowSum(sparse.counts), times=1)
                arma.time[it] <- timeExprs(ArmaRowSum(sparse.counts), times=1)
            }
        }

        writeToFile(Type="sparse (cached)", Ngenes=ngenes, Ncells=ncells, Density=density, 
                    timings=row.time, file="timings_sparse_row.txt", overwrite=overwrite)
        overwrite <- FALSE 
        writeToFile(Type="simple", Ngenes=ngenes, Ncells=ncells, Density=density, 
                    timings=def.time, file="timings_sparse_row.txt", overwrite=overwrite)
        writeToFile(Type="sparse (naive)", Ngenes=ngenes, Ncells=ncells, Density=density, 
                    timings=naive.time, file="timings_sparse_row.txt", overwrite=overwrite)

        if (!skip.arma) { 
            writeToFile(Type="RcppArmadillo", Ngenes=ngenes, Ncells=ncells, Density=density, 
                        timings=arma.time, file="timings_sparse_row.txt", overwrite=overwrite)
            writeToFile(Type="RcppEigen", Ngenes=ngenes, Ncells=ncells, Density=density, 
                        timings=eigen.time, file="timings_sparse_row.txt", overwrite=overwrite)
        }
    }
}

###########################
# Ordered but non-consecutive row access.

ncells <- 1000
density <- 0.01
overwrite <- TRUE
for (ngenes in c(10000, 20000, 50000, 100000)) {
    o <- seq_len(ngenes)
    o <- o[order(o %% 5, o)] # ordered so that there are no consecutive runs.

    naive.time <- row.time <- def.time <- numeric(10)
    for (it in seq_len(10)) { 
        sparse.counts <- rsparsematrix(ngenes, ncells, density)
        dense.counts <- as.matrix(sparse.counts)
        row.time[it] <- timeExprs(BeachmatRowSumRandom(sparse.counts, o), times=1)
        def.time[it] <- timeExprs(BeachmatRowSumRandom(dense.counts, o))
        naive.time[it] <- timeExprs(NaiveSparseRowSumRandom(sparse.counts, o), times=1)
    }

    writeToFile(Type="sparse (cached)", Ngenes=ngenes, Ncells=ncells, Density=density, 
                timings=row.time, file="timings_sparse_row_ordered.txt", overwrite=overwrite)
    overwrite <- FALSE 
    writeToFile(Type="simple", Ngenes=ngenes, Ncells=ncells, Density=density, 
                timings=def.time, file="timings_sparse_row_ordered.txt", overwrite=overwrite)
    writeToFile(Type="sparse (naive)", Ngenes=ngenes, Ncells=ncells, Density=density, 
                timings=naive.time, file="timings_sparse_row_ordered.txt", overwrite=overwrite)
}

###########################
# Random row access.

ncells <- 1000
density <- 0.01
overwrite <- TRUE
for (ngenes in c(10000, 20000, 50000, 100000)) {
    naive.time <- row.time <- def.time <- numeric(10)
    for (it in seq_len(10)) { 
        sparse.counts <- rsparsematrix(ngenes, ncells, density)
        dense.counts <- as.matrix(sparse.counts)
        o <- sample(ngenes)

        row.time[it] <- timeExprs(BeachmatRowSumRandom(sparse.counts, o), times=1)
        def.time[it] <- timeExprs(BeachmatRowSumRandom(dense.counts, o))
        naive.time[it] <- timeExprs(NaiveSparseRowSumRandom(sparse.counts, o), times=1)
    }

    writeToFile(Type="sparse (cached)", Ngenes=ngenes, Ncells=ncells, Density=density, 
                timings=row.time, file="timings_sparse_row_random.txt", overwrite=overwrite)
    overwrite <- FALSE 
    writeToFile(Type="simple", Ngenes=ngenes, Ncells=ncells, Density=density, 
                timings=def.time, file="timings_sparse_row_random.txt", overwrite=overwrite)
    writeToFile(Type="sparse (naive)", Ngenes=ngenes, Ncells=ncells, Density=density, 
                timings=naive.time, file="timings_sparse_row_random.txt", overwrite=overwrite)
}

###########################
# Wrapping up

sessionInfo()

