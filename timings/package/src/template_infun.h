#ifndef BEACHTIME_TEMPLATE_INFUN_H
#define BEACHTIME_TEMPLATE_INFUN_H

/* Calculate the row and column sums after access. */

template <class T, class M>  // M is automatically deduced.
Rcpp::NumericVector get_margins(M ptr, const Rcpp::IntegerVector& mode) {
    if (mode.size()!=1) { 
        throw std::runtime_error("'mode' should be an integer scalar"); 
    }
    const int Mode=mode[0];
    const size_t& nrows=ptr->get_nrow();
    const size_t& ncols=ptr->get_ncol();

    if (Mode==1) { 
        // By column.
        Rcpp::NumericVector output(ncols);
        T target(nrows);
        for (int c=0; c<ncols; ++c) {
            ptr->get_col(c, target.begin());
            output[c]=std::accumulate(target.begin(), target.end(), 0.0);
        }
        return output;
    } else if (Mode==2) { 
        // By row.
        Rcpp::NumericVector output(nrows);
        T target(ncols);
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

/* Test random row/column access. */

template <class T, class M>  
Rcpp::NumericVector get_random_margins(M ptr, const Rcpp::IntegerVector& mode, const Rcpp::IntegerVector& order) {
    if (mode.size()!=1) { 
        throw std::runtime_error("'mode' should be an integer scalar"); 
    }
    const int Mode=mode[0];
    const size_t& nrows=ptr->get_nrow();
    const size_t& ncols=ptr->get_ncol();

    if (Mode==1) { 
        // By column.
        Rcpp::NumericVector output(ncols);
        if (order.size()!=ncols) {
            throw std::runtime_error("order vector should have length equal to the number of columns");
        }
        T target(nrows);

        for (auto oIt=order.begin(); oIt!=order.end(); ++oIt) {
            ptr->get_col(*oIt, target.begin()); // assume zero-indexed.
            output[*oIt]=std::accumulate(target.begin(), target.end(), 0.0);
        }
        return output;
    } else if (Mode==2) { 
        // By row.
        Rcpp::NumericVector output(nrows);
        if (order.size()!=nrows) {
            throw std::runtime_error("order vector should have length equal to the number of rows");
        }
        T target(ncols);

        for (auto oIt=order.begin(); oIt!=order.end(); ++oIt) {
            ptr->get_row(*oIt, target.begin());
            output[*oIt]=std::accumulate(target.begin(), target.end(), 0.0);
        }
        return output;
    } else { 
        throw std::runtime_error("'mode' should be in [1,2]"); 
    }
}

#endif
