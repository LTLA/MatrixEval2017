#ifndef BEACHTIME_TEMPLATE_INFUN_H
#define BEACHTIME_TEMPLATE_INFUN_H

/* Calculate the row and column sums after access. */

template <class M>  // M is automatically deduced.
Rcpp::NumericVector get_margins(M ptr, const Rcpp::IntegerVector& mode) {
    if (mode.size()!=1) { 
        throw std::runtime_error("'mode' should be an integer scalar"); 
    }
    const int Mode=mode[0];
    const size_t& nrows=ptr->get_nrow();
    const size_t& ncols=ptr->get_ncol();

    if (Mode==1) { 
        // By column.
        Rcpp::NumericVector output(ncols), target(nrows);
        for (int c=0; c<ncols; ++c) {
            ptr->get_col(c, target.begin());
            output[c]=std::accumulate(target.begin(), target.end(), 0.0);
        }
        return output;
    } else if (Mode==2) { 
        // By row.
        Rcpp::NumericVector output(nrows), target(ncols);
        for (int r=0; r<nrows; ++r) {
            ptr->get_row(r, target.begin());
            output[r]=std::accumulate(target.begin(), target.end(), 0.0);
        }
        return output;
    } else { 
        throw std::runtime_error("'mode' should be in [1,2]"); 
    }
}

/* Uses almost the same code, replacing beachmat with Rcpp. */

template <class M>  
Rcpp::NumericVector get_default_margins(M mat, const Rcpp::IntegerVector& mode) {
    if (mode.size()!=1) { 
        throw std::runtime_error("'mode' should be an integer scalar"); 
    }
    const int Mode=mode[0];
    const int nrows=mat.nrow();
    const int ncols=mat.ncol();

    if (Mode==1) { 
        // By column.
        Rcpp::NumericVector output(ncols);
        for (int c=0; c<ncols; ++c) {
            auto target=mat.column(c);
            output[c]=std::accumulate(target.begin(), target.end(), 0.0);
        }
        return output;
    } else if (Mode==2) { 
        // By row.
        Rcpp::NumericVector output(nrows);
        for (int r=0; r<nrows; ++r) {
            auto target=mat.row(r);
            output[r]=std::accumulate(target.begin(), target.end(), 0.0);
        }
        return output;
    } else { 
        throw std::runtime_error("'mode' should be in [1,2]"); 
    }
}

/* Calculate the column sums, constant-style. */

template <class M>  // M is automatically deduced.
Rcpp::NumericVector get_col_margins_const(M ptr) {
    const size_t& nrows=ptr->get_nrow();
    const size_t& ncols=ptr->get_ncol();

    Rcpp::NumericVector output(ncols), target(nrows);
    for (int c=0; c<ncols; ++c) {
        auto out=ptr->get_const_col_indexed(c, target.begin());
        auto n=std::get<0>(out);
        auto it=std::get<2>(out);
        output[c]=std::accumulate(it, it+n, 0.0);
    }
    return output;
}

/* Calculate the row and column sums after RANDOM access. Arguably we could do this in get_margins(), 
 * but I don't want to have to modify the functions to accept a 0:(n-1) ordering vector when they don't have to.
 */

template <class M>  // M is automatically deduced.
Rcpp::NumericVector get_margins_random(M ptr, const Rcpp::IntegerVector& mode, const Rcpp::IntegerVector order) {
    if (mode.size()!=1) { 
        throw std::runtime_error("'mode' should be an integer scalar"); 
    }
    const int Mode=mode[0];
    const size_t& nrows=ptr->get_nrow();
    const size_t& ncols=ptr->get_ncol();

    if (Mode==1) { 
        // By column, assuming 'order' is in [0, ncols).
        Rcpp::NumericVector output(order.size()), target(nrows);
        auto oIt=output.begin();

        for (const auto& c : order) {
            ptr->get_col(c, target.begin());
            (*oIt)=std::accumulate(target.begin(), target.end(), 0.0);
            ++oIt;
        }
        return output;
    } else if (Mode==2) { 
        // By row, assuming 'order' is in [0, nrows).
        Rcpp::NumericVector output(order.size()), target(ncols);
        auto oIt=output.begin();

        for (const auto& r : order) {
            ptr->get_row(r, target.begin());
            (*oIt)=std::accumulate(target.begin(), target.end(), 0.0);
            ++oIt;
        }
        return output;
    } else { 
        throw std::runtime_error("'mode' should be in [1,2]"); 
    }
}


#endif
