#ifndef BEACHTEST_H
#define BEACHTEST_H

#include "Rcpp.h"
#include "beachmat/numeric_matrix.h"

extern "C" { 

SEXP get_numeric_margins(SEXP, SEXP);

SEXP get_numeric_default_margins(SEXP, SEXP);

SEXP get_numeric_random_margins(SEXP, SEXP, SEXP);

//SEXP get_simple_margins(SEXP);

}

#endif
