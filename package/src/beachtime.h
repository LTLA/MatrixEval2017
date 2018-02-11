#ifndef BEACHTEST_H
#define BEACHTEST_H

#include <RcppArmadillo.h>
#include <RcppEigen.h>
#include "beachmat/numeric_matrix.h"

extern "C" { 

SEXP get_numeric_margins(SEXP, SEXP);

SEXP get_numeric_default_margins(SEXP, SEXP);

SEXP get_numeric_sparse_row_margins(SEXP);

SEXP get_numeric_col_margins_const(SEXP);

SEXP get_numeric_simple_row_margins(SEXP);

SEXP get_numeric_margins_random(SEXP, SEXP, SEXP);

SEXP get_numeric_sparse_row_margins_random(SEXP, SEXP);

SEXP get_numeric_arma_margins(SEXP, SEXP);

SEXP get_numeric_eigen_margins(SEXP, SEXP);

SEXP edit_numeric_matrix(SEXP);

SEXP edit_numeric_arma_matrix(SEXP);

}

#endif
