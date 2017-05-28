# Simulating row and column access with a base matrix.

library(beachtime)

###########################
# Consecutive column access

ngenes <- 10000    
overwrite <- TRUE
for (ncells in c(1000, 2000, 5000, 10000)) {

    col.time <- def.time <- numeric(10)
    for (it in seq_len(10)) {         
        dense.counts <- matrix(rnorm(ngenes*ncells), ngenes, ncells)
        col.time[it] <- timeExprs(BeachmatColSum(dense.counts))
        def.time[it] <- timeExprs(DefaultColSum(dense.counts))
    }

    write.table(data.frame(Access="col", Type="beachmat", Ngenes=ngenes, Ncells=ncells, Time=mean(col.time), SE=se(col.time)),
                file="timings_base_col.txt", append=!overwrite, col.names=overwrite, sep="\t", quote=FALSE, row.names=FALSE)
    overwrite <- FALSE 
    write.table(data.frame(Access="col", Type="Rcpp", Ngenes=ngenes, Ncells=ncells, Time=mean(def.time), SE=se(def.time)),
                file="timings_base_col.txt", append=!overwrite, col.names=overwrite, sep="\t", quote=FALSE, row.names=FALSE)
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

    write.table(data.frame(Access="row", Type="beachmat", Ngenes=ngenes, Ncells=ncells, Time=mean(row.time), SE=se(row.time)),
                file="timings_base_row.txt", append=!overwrite, col.names=overwrite, sep="\t", quote=FALSE, row.names=FALSE)
    overwrite <- FALSE 
    write.table(data.frame(Access="row", Type="Rcpp", Ngenes=ngenes, Ncells=ncells, Time=mean(def.time), SE=se(def.time)),
                file="timings_base_row.txt", append=!overwrite, col.names=overwrite, sep="\t", quote=FALSE, row.names=FALSE)
}

