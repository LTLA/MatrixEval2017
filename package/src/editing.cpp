#include "beachtime.h"

SEXP edit_numeric_matrix(SEXP in) {
    Rcpp::NumericMatrix mat(in);
    if (mat.nrow()<2 || mat.ncol()<2) {
        throw std::runtime_error("must be at least a 2x2 matrix");
    }
    auto mrow=mat.row(0);
    mrow[0]=123456;
    auto mcol=mat.column(1);
    mcol[1]=987654;
    return mat;
}

SEXP edit_numeric_arma_matrix(SEXP M) {
    arma::sp_mat res=Rcpp::as<arma::sp_mat>(Rcpp::RObject(M));
    if (res.n_rows<2 || res.n_cols<2) {
        throw std::runtime_error("must be at least a 2x2 matrix");
    }
    auto c=res.col(0);
    c[0]=123456;
    auto r=res.col(1);
    r[1]=987654;
    res.print();
    return Rf_ScalarLogical(1);     
}

