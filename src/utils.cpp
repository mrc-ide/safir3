/* --------------------------------------------------------------------------------
 *  utilities for safir
 *  Sean L. Wu (slwood89@gmail.com)
 *  June 2021
 -------------------------------------------------------------------------------- */

#include "../inst/include/utils.hpp"

// [[Rcpp::export]]
bool compare_floats(const double a, const double b) {
  static double eps = std::sqrt(std::numeric_limits<double>::epsilon());
  return std::fabs(a - b) < eps;
};


//' @title Cross tabulate two vectors with given margins
//' @description this is a replacement for \code{\link[base]{table}} that allows empty
//' cells because the margins have been specified. The input vectors \code{a} and \code{b}
//' must be the same length, this function does no argument checking.
//' @param a one set of observations
//' @param b another set of observations
//' @param a_margin number of distinct values of a (rows)
//' @param b_margin number of distinct values of b (cols)
//' @examples
//' a <- 1:5
//' b <- c(1,2,3,1,2)
//' cross_tab_margins(a,b,5,3)
//' table(a,b)
//' @export
// [[Rcpp::export]]
Rcpp::IntegerMatrix cross_tab_margins(
    const Rcpp::IntegerVector& a,
    const Rcpp::IntegerVector& b,
    const int a_margin,
    const int b_margin
) {
  Rcpp::IntegerMatrix out(a_margin, b_margin);

  for (auto i = 0u; i < a.size(); ++i) {
    out(a[i]-1, b[i]-1)++;
  }

  return out;
};


//' @title Cross tabulate doses and age
//' @description The input vectors \code{doses} and \code{age}
//' must have the same number of values, this function does no argument checking.
//' @param doses a \code{\link[individual]{IntegerVariable}}
//' @param age a \code{\link[individual]{IntegerVariable}}
//' @param num_doses number of doses
//' @param num_ages number of age groups
//' @examples
//' \dontrun{
//' a <- IntegerVariable$new(0:4)
//' b <- IntegerVariable$new(c(1,2,3,1,2))
//' cross_tab_doses_age(a$.variable,b$.variable,4,3)
//' table(a$get_values(), b$get_values())
//' }
//' @export
// [[Rcpp::export]]
Rcpp::IntegerMatrix cross_tab_doses_age(
    Rcpp::XPtr<IntegerVariable> doses,
    Rcpp::XPtr<IntegerVariable> age,
    const size_t num_doses,
    const size_t num_ages
) {
  Rcpp::IntegerMatrix out(num_doses + 1, num_ages);

  for (auto i = 0u; i < doses->size; ++i) {
    out(doses->values[i], age->values[i]-1)++;
  }

  return out;
};

//' @title Cross tabulate compartments and age
//' @param compartments a [individual::CategoricalVariable]
//' @param age a [individual::IntegerVariable]
//' @param num_ages number of age groups
//' @param compartment_names a vector giving category names of the [individual::CategoricalVariable]
// [[Rcpp::export]]
Rcpp::IntegerMatrix cross_tab_compartments_age(
    Rcpp::XPtr<CategoricalVariable> compartments,
    Rcpp::XPtr<IntegerVariable> age,
    const size_t num_ages,
    const std::vector<std::string>& compartment_names
) {
  size_t ncompartment{compartment_names.size()};
  Rcpp::IntegerMatrix out(num_ages, ncompartment);

  for (auto j = 0u; j < ncompartment; ++j) {
    individual_index_t comp_bset = compartments->get_index_of(compartment_names[j]);
    for (auto i = 1; i <= num_ages; ++i) {
      individual_index_t age_bset = age->get_index_of_set(i);
      age_bset &= comp_bset;
      out(i - 1, j) = age_bset.size();
    }
  }

  return out;
};


// [[Rcpp::export]]
Rcpp::NumericMatrix cross_tab_margins_internal(
    const std::vector<int>& a,
    const std::vector<int>& b,
    const int a_margin,
    const int b_margin
) {
  Rcpp::NumericMatrix out(a_margin, b_margin);

  for (auto i = 0u; i < a.size(); ++i) {
    out(a[i]-1, b[i]-1) += 1.0;
  }

  return out;
};

//' @title Tabulate a vector of observations
//' @description Tabulate a vector \code{a} whose values fall into a set of integers
//' of maximum value \code{nbins}. This function does no argument checking so please
//' ensure the maximum value of observations is not greater than \code{nbins}.
//' @param a a set of observations
//' @param nbins number of bins
//' @examples
//' nbin <- 10
//' a <- sample.int(n = nbin,size = 100,replace = TRUE)
//' tabulate(bin = a,nbins = nbin)
//' tab_bins(a = a,nbins = nbin)
//' @export
// [[Rcpp::export]]
std::vector<double> tab_bins(
    const std::vector<int>& a,
    const int nbins
) {
  std::vector<double> out(nbins, 0);

  for (auto i = 0u; i < a.size(); ++i) {
    out[a[i]-1] += 1.0;
  }

  return out;
};

//' @title Tabulate a weighted vector of observations
//' @description Similar to [safir::tab_bins] but instead of each observation
//' being implicitly given weight 1, it now has weight given by `wt[i]`.
//' @param a a set of observations
//' @param wt a set of weights
//' @param nbins number of bins
//' @export
// [[Rcpp::export]]
std::vector<double> tab_bins_weighted(
    const std::vector<int>& a,
    const std::vector<double>& wt,
    const int nbins
) {
  std::vector<double> out(nbins, 0);
  if (a.size() != wt.size()) {
    Rcpp::stop("vectors 'a' and 'wt' must be the same length!");
  }

  for (auto i = 0u; i < a.size(); ++i) {
    out[a[i]-1] += wt[i];
  }

  return out;
};

//' @title Get contact matrix
//' @description Get the contact matrix at some specific day (1st dimension of array).
//' @param array the mixing matrix array (days x age x age)
//' @param i the day (indexes the first dimension, assumes zero indexing)
//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix get_contact_matrix_cpp(
    SEXP array
) {
  SEXP dims = Rf_getAttrib(array, R_DimSymbol);
  int d1 = INTEGER(dims)[0];
  int d2 = INTEGER(dims)[1];

  Rcpp::NumericMatrix out(d1,d2);

  double* array_ptr = REAL(array);

  for (auto i = 0u; i < d1; ++i) {
    for (auto j = 0u; j < d2; ++j) {
      out(i,j) = array_ptr[(i * d2) + j];
    }
  }

  return out;
};

//' @title Get a value from a vector
//' @description Get a value at some specific day
//' @param vector_set the set of values
//' @param i the day (assumes zero indexing)
//' @export
// [[Rcpp::export]]
double get_vector_cpp(
    SEXP vector_set,
    const size_t i
){
  double* beta_set_ptr = REAL(vector_set);
  return beta_set_ptr[i];
};

//' @title Multiply a matrix by a integer vector
//' @param m a matrix
//' @param a a vector (must have length equal to number of columns of \code{m})
//' @export
// [[Rcpp::export]]
std::vector<double> matrix_vec_mult_cpp(
    const Rcpp::NumericMatrix& m,
    const std::vector<double>& a
) {
  std::vector<double> out(a.size(), 0.);

  for (auto j = 0u; j < m.ncol(); ++j) {
    for (auto i = 0u; i < m.ncol(); ++i) {
      out[i] += m(i,j) * a[j];
    }
  }

  return out;
};

//' @title Multiply a matrix by a integer vector and a double vector
//' @param m a matrix
//' @param a a vector of double (must have length equal to number of columns of \code{m})
//' @param b a vector of double (must have length equal to number of columns of \code{m})
//' @export
// [[Rcpp::export]]
std::vector<double> matrix_2vec_mult_cpp(
    const Rcpp::NumericMatrix& m,
    const std::vector<double>& a,
    const std::vector<double>& b
) {
  std::vector<double> out(a.size(), 0.);

  for (auto j = 0u; j < m.ncol(); ++j) {
    for (auto i = 0u; i < m.ncol(); ++i) {
      out[i] += m(i,j) * a[j] * b[j];
    }
  }

  return out;
};

//' @title Element-wise multiply two matrices and take row sums
//' @description This function does no argument checking, please make sure \code{a} and \code{b}
//' are matrices of the same dimension.
//' @param a a matrix
//' @param b a matrix
//' @export
// [[Rcpp::export]]
std::vector<double> mult_2matrix_rowsum(
    const Rcpp::NumericMatrix& a,
    const Rcpp::NumericMatrix& b
) {
  std::vector<double> out(a.nrow(), 0.);
  for (auto i = 0u; i < a.nrow(); ++i) {
    for (auto j = 0u; j < a.ncol(); ++j) {
      out[i] += a(i,j) * b(i,j);
    }
  }
  return out;
};

// [[Rcpp::export]]
double get_proportion_vaccinated_nimue_internal(
    Rcpp::XPtr<IntegerVariable> discrete_age,
    Rcpp::XPtr<individual_index_t> vaccinated,
    const int age
) {
  individual_index_t age_bset = discrete_age->get_index_of_set(age);
  double N = age_bset.size();
  individual_index_t vaccinated_bset(*vaccinated);
  vaccinated_bset &= age_bset;
  return static_cast<double>(vaccinated_bset.size()) / N;
};
