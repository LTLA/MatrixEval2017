# Functions for summation of random-access rows or columns.

BeachmatRowSumRandom <- function(M, order) {
    .Call(cxx_get_numeric_margins_random, M, 2L, order - 1L)
}

BeachmatColSumRandom <- function(M, order) {
    .Call(cxx_get_numeric_margins_random, M, 1L, order - 1L)
}

NaiveSparseRowSumRandom <- function(M, order) {
    .Call(cxx_get_numeric_sparse_row_margins_random, M, order - 1L)
}

