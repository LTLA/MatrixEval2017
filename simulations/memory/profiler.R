library(Matrix)
library(HDF5Array)
overwrite <- TRUE
ngenes <- 10000
ncells <- 1000
dir.create("memprof")

# For dense counts:

library(beachtime)
dense.counts <- matrix(rnorm(ngenes*ncells), ngenes, ncells)
gc()
Rprofmem("memprof/dense_col.prof")
x <- BeachmatColSum(dense.counts)
Rprofmem(NULL)

gc()
Rprofmem("memprof/dense_row.prof")
x <- BeachmatRowSum(dense.counts)
Rprofmem(NULL)

# For HDF5 counts:

hdf5.counts <- as(dense.counts, "HDF5Matrix")
gc()
Rprofmem("memprof/hdf5_col.prof")
x <- BeachmatColSum(hdf5.counts)
Rprofmem(NULL)

gc()
Rprofmem("memprof/hdf5_row.prof")
x <- BeachmatRowSum(hdf5.counts)
Rprofmem(NULL)

# For sparse counts:

sparse.counts <- rsparsematrix(ngenes, ncells, 0.01)
gc()
Rprofmem("memprof/sparse_col.prof")
x <- BeachmatColSum(sparse.counts)
Rprofmem(NULL)

gc()
Rprofmem("memprof/sparse_row.prof")
x <- BeachmatRowSum(sparse.counts)
Rprofmem(NULL)



