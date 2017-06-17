a <- matrix(0, 2, 2)
.Call(beachtime:::cxx_edit_numeric_matrix, a)

library(Matrix)
a <- as(diag(2), "dgCMatrix")
.Call(beachtime:::cxx_edit_numeric_arma_matrix, a)
