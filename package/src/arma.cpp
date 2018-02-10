#include "beachtime.h"

SEXP get_numeric_arma_margins (SEXP M, SEXP mode) {
    BEGIN_RCPP
    arma::sp_mat res=Rcpp::as<arma::sp_mat>(Rcpp::RObject(M));

    // Getting row or column sums.
    Rcpp::IntegerVector stuff(mode);
    if (stuff.size()!=1) { 
        throw std::runtime_error("'mode' should be an integer scalar"); 
    }
    const int Mode=stuff[0];
    const int nrows=res.n_rows;
    const int ncols=res.n_cols;

    if (Mode==1) { 
        // By column.
        Rcpp::NumericVector output(ncols);
        for (int c=0; c<ncols; ++c) {
            auto target=res.col(c);
            output[c]=std::accumulate(target.begin(), target.end(), 0.0);
        }
        return output;
    } else if (Mode==2) { 
        // By row.
        Rcpp::NumericVector output(nrows);
        for (int r=0; r<nrows; ++r) {
            auto target=res.row(r);
            output[r]=std::accumulate(target.begin(), target.end(), 0.0);
        }
        return output;
    } else { 
        throw std::runtime_error("'mode' should be in [1,2]"); 
    }
    END_RCPP
}
