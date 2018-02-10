library(pkgDepTools)
library(BiocInstaller)
biocUrl <- biocinstallRepos()["BioCsoft"]
biocDeps <- makeDepGraph(biocUrl, type="source", dosize=FALSE)
biocDepsOnMe <- reverseEdgeDirections(biocDeps)

# Uses Rcpp.
usesRcpp <- dijkstra.sp(biocDepsOnMe, start="Rcpp")$distance
usesRcpp <- usesRcpp[is.finite(usesRcpp)]
length(usesRcpp) - 1 ## don

usesRcppEigen <- dijkstra.sp(biocDepsOnMe, start="RcppEigen")$distance
usesRcppEigen <- usesRcppEigen[is.finite(usesRcppEigen)]
length(usesRcppEigen) - 1 ## don

usesRcppArmadillo <- dijkstra.sp(biocDepsOnMe, start="RcppArmadillo")$distance
usesRcppArmadillo <- usesRcppArmadillo[is.finite(usesRcppArmadillo)]
length(usesRcppArmadillo) - 1 ## don
