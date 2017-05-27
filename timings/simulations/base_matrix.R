# Simulating row and column access with a base matrix.

library(beachtime)
library(microbenchmark)
se <- function(x) { sd(x)/sqrt(length(x)) }

ngenes <- 10000    
overwrite <- TRUE
for (ncells in c(1000, 2000, 5000, 10000)) {
    dense.counts <- matrix(rnorm(ngenes*ncells), ngenes, ncells)
    col.time <- microbenchmark(BeachmatColSum(dense.counts), times=20)$time/1e6 # milliseconds.
    def.time <- microbenchmark(DefaultColSum(dense.counts), times=20)$time/1e6

    write.table(data.frame(Access="col", Type="beachmat", Ngenes=ngenes, Ncells=ncells, Time=mean(col.time), SE=se(col.time)),
                file="timings_base_col.txt", append=!overwrite, col.names=overwrite, sep="\t", quote=FALSE, row.names=FALSE)
    overwrite <- FALSE 
    write.table(data.frame(Access="def", Type="Rcpp", Ngenes=ngenes, Ncells=ncells, Time=mean(def.time), SE=se(def.time)),
                file="timings_base_col.txt", append=!overwrite, col.names=overwrite, sep="\t", quote=FALSE, row.names=FALSE)
}

ncells <- 1000
overwrite <- TRUE
for (ngenes in c(10000, 20000, 50000, 100000)) {
    dense.counts <- matrix(rnorm(ngenes*ncells), ngenes, ncells)
    row.time <- microbenchmark(BeachmatRowSum(dense.counts), times=20)$time/1e6 # milliseconds.
    def.time <- microbenchmark(DefaultRowSum(dense.counts), times=20)$time/1e6

    write.table(data.frame(Access="row", Type="beachmat", Ngenes=ngenes, Ncells=ncells, Time=mean(row.time), SE=se(row.time)),
                file="timings_base_row.txt", append=!overwrite, col.names=overwrite, sep="\t", quote=FALSE, row.names=FALSE)
    overwrite <- FALSE 
    write.table(data.frame(Access="def", Type="Rcpp", Ngenes=ngenes, Ncells=ncells, Time=mean(def.time), SE=se(def.time)),
                file="timings_base_row.txt", append=!overwrite, col.names=overwrite, sep="\t", quote=FALSE, row.names=FALSE)
}

