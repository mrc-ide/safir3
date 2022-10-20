/* --------------------------------------------------------------------------------
 *  efficacy functions for vaccination
 *  Sean L. Wu (slwood89@gmail.com)
 *  August 2021
 -------------------------------------------------------------------------------- */

#ifndef UTILS_EFFICACY_HPP
#define UTILS_EFFICACY_HPP

#include <limits>

#include <Rcpp.h>
#include <individual.h>

std::vector<double> nat_efficacy_infection_cpp(
    const std::vector<double>& ab_titre,
    const Rcpp::List& parameters,
    const size_t day
);

std::vector<double> nat_efficacy_severe_cpp(
        const std::vector<double>& ab_titre,
        const std::vector<double>& ef_infection,
        const Rcpp::List& parameters,
        const size_t day
);

std::vector<double> nat_efficacy_transmission_cpp(
        const std::vector<double>& ab_titre,
        const Rcpp::List& parameters,
        const size_t day
);

#endif
