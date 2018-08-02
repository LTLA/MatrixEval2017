# Rcpp should read by reference with no copies of the matrix.
# Modifying the accessed values will modify the original matrix.
a <- matrix(0, 2, 2)
.Call(beachtime:::cxx_edit_numeric_matrix, a)
a

# This is also true for Matrix->Vector conversions in Rcpp.
# These share memory, so modifying the Vector object will modify the original matrix.
b <- matrix(0, 2, 2)
.Call(beachtime:::cxx_edit_numeric_vector, b)
b

# RcppArmadillo will copy the original matrix, avoiding effects from C++ manipulations.
# However, its accessors still read by reference and can be modified.
library(Matrix)
C <- as(diag(2), "dgCMatrix")
.Call(beachtime:::cxx_edit_numeric_arma_matrix, C)
C
