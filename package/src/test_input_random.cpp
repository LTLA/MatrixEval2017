#include "beachtime.h"
#include "template_infun.h"

/* Realized access functions. */

SEXP get_numeric_margins_random (SEXP in, SEXP mode, SEXP order) {
    BEGIN_RCPP
    auto ptr=beachmat::create_numeric_matrix(in);
    return get_margins_random(ptr.get(), mode, order);
    END_RCPP
}

/* Compute sparse row margins with a naive binary search. */

SEXP get_numeric_sparse_row_margins_random (SEXP in, SEXP order) {
    BEGIN_RCPP
    Rcpp::S4 mat(in);
    Rcpp::NumericVector x(mat.slot("x"));
    Rcpp::IntegerVector i(mat.slot("i"));
    Rcpp::IntegerVector p(mat.slot("p"));
    Rcpp::IntegerVector ordering(order);
    
    Rcpp::IntegerVector dims(mat.slot("Dim"));
    const int nrow=dims[0], ncol=dims[1];
    Rcpp::NumericVector output(ordering.size());
    auto oIt=output.begin();

    for (const auto& r : ordering) {
        double& current_out=(*oIt);
        ++oIt;

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
