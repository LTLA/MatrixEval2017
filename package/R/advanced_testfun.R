standardMatrixMultiply <- function(left, right) {
    stopifnot(identical(ncol(left), nrow(right)))
    .Call(cxx_standard_matrix_multiplication, left, right)
}

indexedMatrixMultiply <- function(left, right) {
    stopifnot(identical(ncol(left), nrow(right)))
    .Call(cxx_indexed_matrix_multiplication, left, right)
}

cellularDetectionRate <- function(M) {
    .Call(cxx_get_ngenes_expressed, M)
}

numCellsExpressing <- function(M) {
    .Call(cxx_get_ncells_expressing, M)
}
