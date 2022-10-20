/* --------------------------------------------------------------------------------
 *  utilities for safir
 *  Sean L. Wu (slwood89@gmail.com)
 *  June 2021
 -------------------------------------------------------------------------------- */

#ifndef UTILS_HPP
#define UTILS_HPP

#include <vector>
#include <numeric>      // std::iota
#include <algorithm>    // std::sort, std::stable_sort
#include <limits>
#include <functional>
#include <Rcpp.h>
#include <individual.h>

bool compare_floats(const double a, const double b);

Rcpp::IntegerMatrix cross_tab_margins(
    const Rcpp::IntegerVector& a,
    const Rcpp::IntegerVector& b,
    const int a_margin,
    const int b_margin
);

Rcpp::IntegerMatrix cross_tab_doses_age(
        Rcpp::XPtr<IntegerVariable> doses,
        Rcpp::XPtr<IntegerVariable> age,
        const size_t num_doses,
        const size_t num_ages
);

Rcpp::NumericMatrix cross_tab_margins_internal(
        const std::vector<int>& a,
        const std::vector<int>& b,
        const int a_margin,
        const int b_margin
);

std::vector<double> tab_bins(
    const std::vector<int>& a,
    const int nbins
);

std::vector<double> tab_bins_weighted(
        const std::vector<int>& a,
        const std::vector<double>& wt,
        const int nbins
);

Rcpp::NumericMatrix get_contact_matrix_cpp(
    SEXP array
);

double get_vector_cpp(
    SEXP vector_set,
    const size_t i
);

std::vector<double> matrix_vec_mult_cpp(
    const Rcpp::NumericMatrix& m,
    const std::vector<double>& a
);

std::vector<double> matrix_2vec_mult_cpp(
    const Rcpp::NumericMatrix& m,
    const std::vector<double>& a,
    const std::vector<double>& b
);

std::vector<double> mult_2matrix_rowsum(
    const Rcpp::NumericMatrix& a,
    const Rcpp::NumericMatrix& b
);

double get_proportion_vaccinated_nimue_internal(
    Rcpp::XPtr<IntegerVariable> discrete_age,
    Rcpp::XPtr<individual_index_t> vaccinated,
    const int age
);

#endif
