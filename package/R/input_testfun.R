# Creating functions to check each set of matrix inputs.

BeachmatRowSum <- function(M) {
    .Call(cxx_get_numeric_margins, M, 2L)
}

BeachmatColSum <- function(M) {
    .Call(cxx_get_numeric_margins, M, 1L)
}

BeachmatColSumNoCopy <- function(M) {
    .Call(cxx_get_numeric_col_margins_const, M)
}

DefaultRowSum <- function(M) {
    .Call(cxx_get_numeric_default_margins, M, 2L)
}

DefaultColSum <- function(M) {
    .Call(cxx_get_numeric_default_margins, M, 1L)
}

NaiveSparseRowSum <- function(M) {
    .Call(cxx_get_numeric_sparse_row_margins, M)
}

DirectSimpleRowSum <- function(M) {
    .Call(cxx_get_numeric_simple_row_margins, M)
}

ArmaRowSum <- function(M) {
    .Call(cxx_get_numeric_arma_margins, M, 2L)
}

ArmaColSum <- function(M) {
    .Call(cxx_get_numeric_arma_margins, M, 1L)
}

EigenRowSum <- function(M) {
    .Call(cxx_get_numeric_arma_margins, M, 2L)
}

EigenColSum <- function(M) {
    .Call(cxx_get_numeric_arma_margins, M, 1L)
}
