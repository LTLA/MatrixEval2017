#include "beachtime.h"
#include "template_infun.h"

/* Const column functions. */

SEXP get_numeric_col_margins_const (SEXP in) {
    BEGIN_RCPP
    auto ptr=beachmat::create_numeric_matrix(in);
    return get_col_margins_const(ptr.get());
    END_RCPP
}

/* Compute simple row margins without copying. */

SEXP get_numeric_simple_row_margins (SEXP in) {
    BEGIN_RCPP
    Rcpp::NumericMatrix mat(in);
    const int nrows=mat.nrow();
    const int ncols=mat.ncol();

    Rcpp::NumericVector output(mat.nrow());
    for (int r=0; r<nrows; ++r) {
        auto mIt=mat.begin();
        double& curout=output[r];
        for (int c=0; c<ncols; ++c, mIt += nrows) {
            curout += *mIt;
        }
    }
    
    return output;
    END_RCPP
}


