/* --------------------------------------------------------------------------------
 *  infection process for individual infectiousness transmission model
 -------------------------------------------------------------------------------- */

#include <Rcpp.h>
#include <individual.h>
#include "../inst/include/utils.hpp"

//' @title C++ infection process (modified squire transmission model)
//' @description this is an internal function, you should use the R interface
//' for type checking, \code{\link{infection_process_cpp}}
//' @param parameters a list of parameters from \code{\link{get_parameters}}
//' @param states a \code{\link[individual]{CategoricalVariable}}
//' @param discrete_age a \code{\link[individual]{IntegerVariable}}
//' @param receive_antivirals a \code{\link[individual]{IntegerVariable}}
//' @param infectiousness_start_time a \code{\link[individual]{IntegerVariable}}
//' @param infection a \code{\link[individual]{TargetedEvent}}
//' @param dt size of time step
// [[Rcpp::export]]
Rcpp::XPtr<process_t> infection_process_cpp_internal(
    Rcpp::List parameters,
    Rcpp::XPtr<CategoricalVariable> states,
    Rcpp::XPtr<IntegerVariable> discrete_age,
    Rcpp::XPtr<IntegerVariable> infectiousness_start_time,
    Rcpp::XPtr<IntegerVariable> wearing_ppe,
    Rcpp::XPtr<TargetedEvent> infection,
    const double dt
) {

  // Required parameters
  int number_of_age_groups =
    Rcpp::as<int>(parameters["number_of_age_groups"]);
  double ppe_transmission_multiplier =
    Rcpp::as<double>(parameters["ppe_transmission_multiplier"]);
  std::vector<double> infectiousness_profile =
    Rcpp::as<std::vector<double>>(parameters["infectiousness_profile"]);
  SEXP contact_matrix = parameters["contact_matrix"];
  SEXP beta_set = parameters["beta_set"];
  SEXP lambda_external_vector = parameters["lambda_external"];

  // Contact matrix
  Rcpp::NumericMatrix m = get_contact_matrix_cpp(contact_matrix);

  // Infection process function
  return Rcpp::XPtr<process_t>(
    new process_t([parameters, states, discrete_age, infectiousness_start_time, wearing_ppe,
                   infection, dt, infectiousness_profile, ppe_transmission_multiplier,
                   m, beta_set, lambda_external_vector, number_of_age_groups](size_t t) mutable {

      // Current day (subtract one for zero-based indexing)
      size_t tnow = std::ceil((double)t * dt) - 1.;

      // Force of infection from contact outside the population
      double lambda_external = get_vector_cpp(lambda_external_vector, tnow);

      // Infectious persons
      individual_index_t infectious = states->get_index_of("I");

      // Susceptible persons
      individual_index_t susceptible = states->get_index_of("S");

      if (susceptible.size() > 0) {

        // Force of infection from external contacts
        std::vector<double> lambda(susceptible.size(), lambda_external);

        // Force of infection from transmission
        if (infectious.size() > 0) {

          // Force of infection by age group
          std::vector<double> lambda_age(number_of_age_groups, 1.0);

          // For the infectious individuals, record their age group, when their period of
          //infectiousness started and whether or not they are wearing PPE
          std::vector<int> ages = discrete_age->get_values(infectious);
          std::vector<int> start_time = infectiousness_start_time->get_values(infectious);
          std::vector<int> wearing_ppe_infectious = wearing_ppe->get_values(infectious);

          // Determine the current level of infectiousness for the infectious individuals
          std::vector<double> infectiousness(ages.size(), 0.0);
          int offset, ppe_multiplier;
          for(int n=0; n<start_time.size(); n++){
            offset = (int)((t - start_time[n]) * dt);
            if(offset >= infectiousness_profile.size()){
              infectiousness[n] = 0.0;
            } else {
              ppe_multiplier = 1 + (wearing_ppe_infectious[n] * (ppe_transmission_multiplier - 1));
              infectiousness[n] = infectiousness_profile[offset] * ppe_multiplier;
            }
          }

          // Calculate force of infection by age group
          double beta = get_vector_cpp(beta_set, tnow);
          for(int a=0; a<number_of_age_groups; a++){
            for(int n=0; n<ages.size(); n++){
              lambda_age[a] *= 1 - (dt * m(a, ages[n] - 1) * beta * infectiousness[n]);
            }
          }
          for(int a=0; a<number_of_age_groups; a++){
            lambda_age[a] = 1 - lambda_age[a];
          }

          // For the susceptible individuals, record their age group and whether or not they are
          // wearing PPE
          std::vector<int> sus_ages = discrete_age->get_values(susceptible);
          std::vector<int> wearing_ppe_susceptible = wearing_ppe->get_values(susceptible);

          // Using the force of infection for their age group, calculate the probability of each
          // susceptible being infected
          for (int n=0; n<sus_ages.size(); n++) {
            ppe_multiplier = 1 + (wearing_ppe_susceptible[n] * (ppe_transmission_multiplier - 1));
            lambda[n] += lambda_age[sus_ages[n] - 1] * ppe_multiplier;
          }

        }

        // Sample the susceptible population using the probabilities in lambda
        bitset_sample_multi_internal(susceptible, lambda.begin(), lambda.end());

        // Queue the infection events
        if (susceptible.size() > 0) {
          infection->schedule(susceptible, 0.0);
        }

      }

    }),
    true
  );
};
