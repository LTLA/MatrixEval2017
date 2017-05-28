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

/* Other access functions */

SEXP get_numeric_arma_margins (SEXP M, SEXP mode) {
    BEGIN_RCPP
    
    Rcpp::S4 mat(M);
    Rcpp::IntegerVector dims = mat.slot("Dim");
    arma::urowvec i = Rcpp::as<arma::urowvec>(mat.slot("i"));
    arma::urowvec p = Rcpp::as<arma::urowvec>(mat.slot("p"));     
    arma::vec x     = Rcpp::as<arma::vec>(mat.slot("x"));

    int nrow = dims[0], ncol = dims[1];
    arma::sp_mat res(nrow, ncol);

    // create space for values, and copy
    arma::access::rw(res.values) = arma::memory::acquire_chunked<double>(x.size() + 1);
    arma::arrayops::copy(arma::access::rwp(res.values), x.begin(), x.size() + 1);

    // create space for row_indices, and copy
    arma::access::rw(res.row_indices) = arma::memory::acquire_chunked<arma::uword>(i.size() + 1);
    arma::arrayops::copy(arma::access::rwp(res.row_indices), i.begin(), i.size() + 1);
    
    // create space for col_ptrs, and copy 
    arma::access::rw(res.col_ptrs) = arma::memory::acquire<arma::uword>(p.size() + 2);
    arma::arrayops::copy(arma::access::rwp(res.col_ptrs), p.begin(), p.size() + 1);

    // important: set the sentinel as well
    arma::access::rwp(res.col_ptrs)[p.size()+1] = std::numeric_limits<arma::uword>::max();
    
    // set the number of non-zero elements
    arma::access::rw(res.n_nonzero) = x.size();

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

SEXP get_numeric_simple_margins (SEXP in) {
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
