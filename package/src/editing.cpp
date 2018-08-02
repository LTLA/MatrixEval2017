#include "beachtime.h"

// Checks that the accessors are read-by-reference.
// This means that modifications occur in-place in 'mat' (and hence, in 'in').
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

// This just checks that the memory is recycled between Rcpp::Matrices and Vectors of the same type.
// That is, a deep copy is _not_ created, and the modifications below manifest in the passed matrix.
SEXP edit_numeric_vector(SEXP in) {
    Rcpp::NumericVector vec(in);
    if (vec.size()<2) {
        throw std::runtime_error("must be at least length 2");
    }

    vec[0]=123456;
    vec[1]=987654;
    return vec;
}

// Checks that accessors are read-by-reference.
// This means that any modifications occur in-place within the sp_mat object.
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

