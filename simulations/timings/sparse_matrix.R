# Simulating row and column access with a base matrix.

library(beachtime)
library(Matrix)

###########################
# Consecutive column access

ngenes <- 10000    
overwrite <- TRUE
for (ncells in c(1000, 2000, 5000, 10000)) {
    for (density in c(0.01, 0.05, 0.1, 0.2)) { 
        
        col.time <- arma.time <- def.time <- numeric(10)
        for (it in seq_len(10)) { 
            sparse.counts <- rsparsematrix(ngenes, ncells, density)
            dense.counts <- as.matrix(sparse.counts)
            col.time[it] <- timeExprs(BeachmatColSum(sparse.counts))
            arma.time[it] <- timeExprs(ArmaColSum(sparse.counts))
            def.time[it] <- timeExprs(BeachmatColSum(dense.counts))
        }
    
        writeToFile(Type="sparse", Ngenes=ngenes, Ncells=ncells, Density=density, 
                    timings=col.time, file="timings_sparse_col.txt", overwrite=overwrite)
        overwrite <- FALSE 
        writeToFile(Type="RcppArmadillo", Ngenes=ngenes, Ncells=ncells, Density=density, 
                    timings=arma.time, file="timings_sparse_col.txt", overwrite=overwrite)
        writeToFile(Type="simple", Ngenes=ngenes, Ncells=ncells, Density=density, 
                    timings=def.time, file="timings_sparse_col.txt", overwrite=overwrite)
    }
}

###########################
# Consecutive row access

ncells <- 1000
overwrite <- TRUE
skip.arma <- FALSE
for (ngenes in c(10000, 20000, 50000, 100000)) {
    for (density in c(0.01, 0.05, 0.1, 0.2)) { 

        naive.time <- row.time <- arma.time <- def.time <- numeric(10)
        for (it in seq_len(10)) { 
            sparse.counts <- rsparsematrix(ngenes, ncells, density)
            dense.counts <- as.matrix(sparse.counts)
            row.time[it] <- timeExprs(BeachmatRowSum(sparse.counts))
            def.time[it] <- timeExprs(BeachmatRowSum(dense.counts))
            naive.time[it] <- timeExprs(NaiveSparseRowSum(sparse.counts), times=1) # slow, so only once.

            # This takes some time, so will only run for the lowest number of genes, and once.
            if (!skip.arma) {
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
        }
    }
            
    skip.arma <- TRUE
}

###########################
# Wrapping up

sessionInfo()

