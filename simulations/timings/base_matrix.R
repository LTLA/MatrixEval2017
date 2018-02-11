# Simulating row and column access with a base matrix.

library(beachtime)

###########################
# Consecutive column access

ngenes <- 10000    
overwrite <- TRUE
for (ncells in c(1000, 2000, 5000, 10000)) {

    col.time <- colnocopy.time <- def.time <- numeric(10)
    for (it in seq_len(10)) {         
        dense.counts <- matrix(rnorm(ngenes*ncells), ngenes, ncells)
        col.time[it] <- timeExprs(BeachmatColSum(dense.counts))
        colnocopy.time[it] <- timeExprs(BeachmatColSumNoCopy(dense.counts))
        def.time[it] <- timeExprs(DefaultColSum(dense.counts))
    }

    writeToFile(Type="beachmat", Ngenes=ngenes, Ncells=ncells,
                timings=col.time, file="timings_base_col.txt", overwrite=overwrite)
    overwrite <- FALSE 
    writeToFile(Type="beachmat (no copy)", Ngenes=ngenes, Ncells=ncells,
                timings=colnocopy.time, file="timings_base_col.txt", overwrite=overwrite)
    writeToFile(Type="Rcpp", Ngenes=ngenes, Ncells=ncells, 
                timings=def.time, file="timings_base_col.txt", overwrite=overwrite)
}

###########################
# Consecutive row access

ncells <- 1000
overwrite <- TRUE
for (ngenes in c(10000, 20000, 50000, 100000)) {

    row.time <- def.time <- numeric(10)
    for (it in seq_len(10)) {         
        dense.counts <- matrix(rnorm(ngenes*ncells), ngenes, ncells)
        row.time[it] <- timeExprs(BeachmatRowSum(dense.counts))
        def.time[it] <- timeExprs(DefaultRowSum(dense.counts))
    }

    writeToFile(Type="beachmat", Ngenes=ngenes, Ncells=ncells, 
                timings=row.time, file="timings_base_row.txt", overwrite=overwrite)
    overwrite <- FALSE 
    writeToFile(Type="Rcpp", Ngenes=ngenes, Ncells=ncells, 
                timings=def.time, file="timings_base_row.txt", overwrite=overwrite)
}

###########################
# Wrapping up

sessionInfo()
