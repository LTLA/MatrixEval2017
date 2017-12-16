#include "beachtime.h"
 
SEXP get_numeric_eigen_margins (SEXP M, SEXP mode) {
    BEGIN_RCPP

    // Getting row or column sums.
    Rcpp::IntegerVector stuff(mode);
    if (stuff.size()!=1) {
        throw std::runtime_error("'mode' should be an integer scalar");
    }
    const int Mode=stuff[0];
   
    // Code obtained from http://gallery.rcpp.org/articles/sparse-iterators/.
    typedef Eigen::MappedSparseMatrix<double> MSpMat;
    const MSpMat X = Rcpp::as<MSpMat>(M);
    const int nrows=X.rows();
    const int ncols=X.cols();

    if (Mode==1) {
        // By column.
        Rcpp::NumericVector output(ncols);
        for (int c=0; c<ncols; ++c){
            double& current=output[c];
            for (int r=0; r<nrows; ++r) { 
                current+=X.coeff(r, c);
            }
        }
        return output;
    } else if (Mode==2) {
        // By row.
        Rcpp::NumericVector output(nrows);
        for (int r=0; r<nrows; ++r) {
            double& current=output[r];
            for (int c=0; c<ncols; ++c){
                current+=X.coeff(r, c);
            }
        }
        return output;
    } else {
        throw std::runtime_error("'mode' should be in [1,2]");
    }
    END_RCPP
}
