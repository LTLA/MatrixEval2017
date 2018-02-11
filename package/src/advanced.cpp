#include "beachtime.h"

SEXP dense_matrix_multiplication(SEXP left, SEXP right) {
    BEGIN_RCPP
    auto lptr=beachmat::create_numeric_matrix(left);
    auto rptr=beachmat::create_numeric_matrix(right);

    const size_t NC_left=lptr->get_ncol(); 
    const size_t NR_left=lptr->get_ncol(); 
    const size_t NC_right=rptr->get_ncol(); 
    const size_t NR_right=rptr->get_ncol(); 

    Rcpp::NumericMatrix output(NR_left, NC_right);
    Rcpp::NumericVector left_out(NC_left), right_out(NR_right);
    
    // Fairly naive approach to matrix multiplication.
    for (size_t row_left=0; row_left<NR_left; ++row_left) {
        lptr->get_row(row_left, left_out.begin());
        auto oIt=output.begin() + row_left;
   
        for (size_t col_right=0; col_right<NC_right; ++col_right) {
            auto right_info=rptr->get_const_col_indexed(col_right, right_out.begin());
            auto rnum=std::get<0>(right_info);
            auto rdex=std::get<1>(right_info);
            auto rval=std::get<2>(right_info);

            double& curval=*(oIt);
            for (size_t i=0; i<rnum; ++i) {
                curval+=left_out[*rdex] * (*rval);
                ++rdex;
                ++rval;
            }
            oIt+=NR_left;
        }
    }

    return output;
    END_RCPP
}

SEXP sparse_matrix_multiplication(SEXP left, SEXP right) {
    BEGIN_RCPP
    auto lptr=beachmat::create_numeric_matrix(left);
    auto rptr=beachmat::create_numeric_matrix(right);

    const size_t NC_left=lptr->get_ncol(); 
    const size_t NR_left=lptr->get_ncol(); 
    const size_t NC_right=rptr->get_ncol(); 
    const size_t NR_right=rptr->get_ncol(); 

    Rcpp::NumericMatrix output(NR_left, NC_right);
    Rcpp::NumericVector left_out(NC_left), right_out(NR_right);
    
    // This uses indexing to avoid dealing with zeroes.
    for (size_t col_right=0; col_right<NC_right; ++col_right) {
        auto right_info=rptr->get_const_col_indexed(col_right, right_out.begin());
        auto rnum=std::get<0>(right_info);
        auto rdex=std::get<1>(right_info);
        auto rval=std::get<2>(right_info);

        for (size_t rcounter=0; rcounter<rnum; ++rcounter, ++rdex, ++rval) {
            const auto& col_left=(*rdex);
            const double& cur_val=(*rval);

            auto left_info=lptr->get_const_col_indexed(col_left, left_out.begin());
            auto lnum=std::get<0>(left_info);
            auto ldex=std::get<1>(left_info);
            auto lval=std::get<2>(left_info);

            auto oIt=output.begin() + col_left * NR_left;
            for (size_t i=0; i<lnum; ++i) {
                *(oIt+*ldex)+=(*lval) * cur_val;
                ++ldex;
                ++lval;
            }
        }
    }

    return output;
    END_RCPP
}
