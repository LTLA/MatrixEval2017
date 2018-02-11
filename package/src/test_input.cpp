#include "beachtime.h"
#include "template_infun.h"

/* Realized access functions. */

SEXP get_numeric_margins (SEXP in, SEXP mode) {
    BEGIN_RCPP
    auto ptr=beachmat::create_numeric_matrix(in);
    return get_margins(ptr.get(), mode);
    END_RCPP
}

SEXP get_numeric_default_margins (SEXP in, SEXP mode) {
    BEGIN_RCPP
    Rcpp::NumericMatrix mat(in);
    return get_default_margins(mat, mode);
    END_RCPP
}

/* Compute sparse row margins with a naive binary search. */

SEXP get_numeric_sparse_row_margins (SEXP in) {
    BEGIN_RCPP
    Rcpp::S4 mat(in);
    Rcpp::NumericVector x(mat.slot("x"));
    Rcpp::IntegerVector i(mat.slot("i"));
    Rcpp::IntegerVector p(mat.slot("p"));
    
    Rcpp::IntegerVector dims(mat.slot("Dim"));
    const int nrow=dims[0], ncol=dims[1];
    Rcpp::NumericVector output(nrow);

    for (int r=0; r<nrow; ++r) {
        double& current_out=output[r];

        for (int c=0; c<ncol; ++c) {
            auto first=i.begin() + p[c];
            auto last=i.begin() + p[c+1];
            auto loc=std::lower_bound(first, last, r);
            if (loc!=last && *loc == r) {
                current_out+=x[loc - i.begin()];
            }
        }
    } 

    return output;
    END_RCPP
}
