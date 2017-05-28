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
        writeToFile(Type="arma", Ngenes=ngenes, Ncells=ncells, Density=density, 
                    timings=arma.time, file="timings_sparse_col.txt", overwrite=overwrite)
        writeToFile(Type="dense", Ngenes=ngenes, Ncells=ncells, Density=density, 
                    timings=def.time, file="timings_sparse_col.txt", overwrite=overwrite)
    }
}

###########################
# Consecutive row access

ncells <- 1000
overwrite <- TRUE
for (ngenes in c(10000, 20000, 50000, 100000)) {
    for (density in c(0.01, 0.05, 0.1, 0.2)) { 

        row.time <- arma.time <- def.time <- numeric(10)
        for (it in seq_len(10)) { 
            sparse.counts <- rsparsematrix(ngenes, ncells, density)
            dense.counts <- as.matrix(sparse.counts)
            row.time[it] <- timeExprs(BeachmatRowSum(sparse.counts))
            arma.time[it] <- timeExprs(ArmaRowSum(sparse.counts), times=1)
            def.time[it] <- timeExprs(BeachmatRowSum(dense.counts))
        }

        writeToFile(Type="sparse", Ngenes=ngenes, Ncells=ncells, Density=density, 
                    timings=row.time, file="timings_sparse_row.txt", overwrite=overwrite)
        overwrite <- FALSE 
        writeToFile(Type="arma", Ngenes=ngenes, Ncells=ncells, Density=density, 
                    timings=arma.time, file="timings_sparse_row.txt", overwrite=overwrite)
        writeToFile(Type="dense", Ngenes=ngenes, Ncells=ncells, Density=density, 
                    timings=def.time, file="timings_sparse_row.txt", overwrite=overwrite)
    }
}

###########################
# Random row access

ncells <- 1000
overwrite <- TRUE
for (ngenes in c(10000, 20000, 50000, 100000)) {
    for (density in c(0.01, 0.05, 0.1, 0.2,)) { 

        row.time <- def.time <- numeric(10)
        for (it in seq_len(10)) { 
            sparse.counts <- rsparsematrix(ngenes, ncells, density)
            o <- sample(ngenes, ngenes)
            i <- seq_len(ngenes)
            row.time[it] <- timeExprs(RandomRowSum(sparse.counts, o))
            def.time[it] <- timeExprs(RandomRowSum(sparse.counts, i))
        }

        writeToFile(Type="random", Ngenes=ngenes, Ncells=ncells, Density=density, 
                    timings=row.time, file="timings_sparse_row_rand.txt", overwrite=overwrite)
        overwrite <- FALSE 
        writeToFile(Type="consecutive", Ngenes=ngenes, Ncells=ncells, Density=density, 
                    timings=def.time, file="timings_sparse_row_rand.txt", overwrite=overwrite)
    }
}

