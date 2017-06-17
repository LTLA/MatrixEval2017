#ifndef BEACHTEST_H
#define BEACHTEST_H

#include <RcppArmadillo.h>
#include "beachmat/numeric_matrix.h"

arma::sp_mat create_sparse_arma(SEXP);

extern "C" { 

SEXP get_numeric_margins(SEXP, SEXP);

SEXP get_numeric_default_margins(SEXP, SEXP);

SEXP get_numeric_simple_row_margins(SEXP);

SEXP get_numeric_sparse_row_margins(SEXP);

SEXP get_numeric_arma_margins(SEXP, SEXP);

SEXP edit_numeric_matrix(SEXP);

SEXP edit_numeric_arma_matrix(SEXP);

}

#endif
