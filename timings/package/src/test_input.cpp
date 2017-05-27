#include "beachtime.h"
#include "template_infun.h"

/* Realized access functions. */

SEXP get_numeric_margins (SEXP in, SEXP mode) {
    BEGIN_RCPP
    auto ptr=beachmat::create_numeric_matrix(in);
    return get_margins<Rcpp::NumericVector>(ptr.get(), mode);
    END_RCPP
}

SEXP get_numeric_default_margins (SEXP in, SEXP mode) {
    BEGIN_RCPP
    Rcpp::NumericMatrix mat(in);
    return get_default_margins(mat, mode);
    END_RCPP
}

SEXP get_numeric_random_margins (SEXP in, SEXP mode, SEXP order) {
    BEGIN_RCPP
    auto ptr=beachmat::create_numeric_matrix(in);
    return get_random_margins<Rcpp::NumericVector>(ptr.get(), mode, order);
    END_RCPP
}

/*
SEXP get_simple_margins (SEXP in) {
    BEGIN_RCPP
    Rcpp::NumericMatrix mat(in);
    const int nrows=mat.nrow();
    const int ncols=mat.ncol();

    Rcpp::NumericVector output(mat.nrow());
//    Rcpp::NumericVector tmp(mat.ncol());
    for (int r=0; r<nrows; ++r) {
//        auto tIt=&(tmp[0]);
        auto mIt=mat.begin();
        double& curout=output[r];
        for (int c=0; c<ncols; ++c, mIt += nrows) {
            curout += *mIt;
        }
    }
    
    return output;
    END_RCPP
}
*/
