#include "beachtime.h"

SEXP get_ncells_expressing (SEXP M) {
    BEGIN_RCPP
    auto ptr=beachmat::create_numeric_matrix(M);
    const size_t NC=ptr->get_ncol();
    const size_t NR=ptr->get_nrow();

    Rcpp::IntegerVector output (NR);
    Rcpp::NumericVector target (NR);

    for (size_t c=0; c<NC; ++c) {
        auto info=ptr->get_const_col_indexed(c, target.begin());
        size_t num=std::get<0>(info);
        auto dex=std::get<1>(info);
        auto val=std::get<2>(info);

        for (size_t i=0; i<num; ++i, ++val, ++dex) {
            if (*val > 0) { ++output[*dex]; }
        }
    }

    return output;
    END_RCPP
}

SEXP get_ngenes_expressed (SEXP M) {
    BEGIN_RCPP
    auto ptr=beachmat::create_numeric_matrix(M);
    const size_t NC=ptr->get_ncol();
    const size_t NR=ptr->get_nrow();

    Rcpp::IntegerVector output (NC);
    Rcpp::NumericVector target (NR);

    for (size_t c=0; c<NC; ++c) {
        auto info=ptr->get_const_col_indexed(c, target.begin());
        size_t num=std::get<0>(info);
        auto val=std::get<2>(info);

        int& count=output[c];
        for (size_t i=0; i<num; ++i, ++val) {
            if (*val > 0) { ++count; }
        }
    }

    return output;
    END_RCPP
}


