# This script checks that beachmat is correctly setting the chunk cache.
# This enables rapid consecutive access from rectangular-chunked HDF5 files.
# The idea is to check that access is competitive with pure row/column chunks.

ngenes <- 10000
ncells <- 10000    
dense.counts <- matrix(rnorm(ngenes*ncells), ngenes, ncells)

library(HDF5Array)
out.bycol <- writeHDF5Array(dense.counts, "bycol.h5", name="yyy", chunkdim=c(ngenes, 1), level=6)
out.byrow <- writeHDF5Array(dense.counts, "byrow.h5", name="yyy", chunkdim=c(1, ncells), level=6)
out.rect <- writeHDF5Array(dense.counts, "rect.h5", name="yyy", chunkdim=c(100, 100), level=6)

out.bycol <- HDF5Array("bycol.h5", "yyy", "double")
out.byrow <- HDF5Array("byrow.h5", "yyy", "double")
out.rect <- HDF5Array("rect.h5", "yyy", "double")

library(beachtime)
timeExprs(beachtime::BeachmatColSum(out.bycol), times=1)
timeExprs(beachtime::BeachmatColSum(out.rect), times=1)
timeExprs(beachtime::BeachmatRowSum(out.byrow), times=1)
timeExprs(beachtime::BeachmatRowSum(out.rect), times=1)

#library(scater)
#scater:::rechunk(out.bycol, 500, byrow=TRUE)
