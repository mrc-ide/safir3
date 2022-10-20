/* --------------------------------------------------------------------------------
 *  infection process for individual infectiousness transmission model
 -------------------------------------------------------------------------------- */

#include <Rcpp.h>
#include <individual.h>
#include "../inst/include/utils.hpp"
#include "../inst/include/efficacy_nat.hpp"

// types for lambda functions used for calculations below
using get_trans_eff_func =
  std::function<std::vector<double>(const individual_index_t&, Rcpp::List, const size_t)>;
using calculate_nat_func =
  std::function<std::vector<double>(const individual_index_t&)>;

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
Rcpp::XPtr<process_t> infection_process_vaccine_cpp_internal(
    Rcpp::List parameters,
    Rcpp::List variables,
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

  // Calculate NAT
  calculate_nat_func calculate_nat;

  Rcpp::Environment ab_titre_R6 = Rcpp::as<Rcpp::Environment>(variables["ab_titre"]);
  Rcpp::XPtr<DoubleVariable> ab_titre(Rcpp::as<SEXP>(ab_titre_R6[".variable"]));

  if (variables.containsElementNamed("ab_titre_inf")) {

    Rcpp::Environment ab_titre_inf_R6 = Rcpp::as<Rcpp::Environment>(variables["ab_titre_inf"]);
    Rcpp::XPtr<DoubleVariable> ab_titre_inf(Rcpp::as<SEXP>(ab_titre_inf_R6[".variable"]));

    calculate_nat =
      [ab_titre, ab_titre_inf](const individual_index_t& index) -> std::vector<double> {

      std::vector<double> nat_vaccine = ab_titre->get_values(index);
      std::vector<double> nat_infection = ab_titre_inf->get_values(index);

      std::vector<double> nat(index.size());

      for (auto i = 0u; i < index.size(); ++i) {
        nat[i] = std::exp(nat_vaccine[i]) + std::exp(nat_infection[i]);
        nat[i] = std::log(nat[i]);
      }

      return nat;
    };
  } else {
    calculate_nat =
      [ab_titre](const individual_index_t& index) -> std::vector<double> {
      return ab_titre->get_values(index);
    };
  }

  // Calculate infectiousness weights
  get_trans_eff_func get_trans_eff;

  if (Rcpp::as<bool>(parameters["nt_efficacy_transmission"])) {
    get_trans_eff = [discrete_age, calculate_nat](const individual_index_t& infectious_bset,
                              Rcpp::List parameters, const size_t timestep) -> std::vector<double> {
      std::vector<int> ages = discrete_age->get_values(infectious_bset);
      std::vector<double> nat_values = calculate_nat(infectious_bset);
      std::vector<double> transmission_efficacy =
        nat_efficacy_transmission_cpp(nat_values, parameters, timestep);
      return(transmission_efficacy);
    };
  } else {
    get_trans_eff = [discrete_age](const individual_index_t& infectious_bset,
                              Rcpp::List parameters, const size_t timestep) -> std::vector<double> {
      std::vector<double> transmission_efficacy(infectious_bset.size(), 1.0);
      return(transmission_efficacy);
    };
  }

  // Infection process function
  return Rcpp::XPtr<process_t>(
    new process_t([parameters, states, discrete_age, infectiousness_start_time, wearing_ppe,
                   infection, dt, infectiousness_profile, ppe_transmission_multiplier,
                   m, beta_set, lambda_external_vector, number_of_age_groups,
                   get_trans_eff, calculate_nat](size_t t) mutable {

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
          // infectiousness started, whether or not they are wearing PPE, and their current
          // reduction in transmission due to previous infection or vaccination
          std::vector<int> ages = discrete_age->get_values(infectious);
          std::vector<int> start_time = infectiousness_start_time->get_values(infectious);
          std::vector<int> wearing_ppe_infectious = wearing_ppe->get_values(infectious);
          std::vector<double> transmission_efficacy = get_trans_eff(infectious, parameters, tnow);

          // Determine the current level of infectiousness for the infectious individuals
          std::vector<double> infectiousness(ages.size(), 0.0);
          int offset, ppe_multiplier;
          for(int n=0; n<start_time.size(); n++){
            offset = (int)((t - start_time[n]) * dt);
            if(offset >= infectiousness_profile.size()){
              infectiousness[n] = 0.0;
            } else {
              ppe_multiplier = 1 + (wearing_ppe_infectious[n] * (ppe_transmission_multiplier - 1));
              infectiousness[n] =
                infectiousness_profile[offset] * ppe_multiplier * transmission_efficacy[n];
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

          // Get NAT efficacy against infection
          std::vector<double> ab_titre_susceptible = calculate_nat(susceptible);
          std::vector<double> infection_efficacy =
            nat_efficacy_infection_cpp(ab_titre_susceptible, parameters, tnow);

          // Using the force of infection for their age group, calculate the probability of each
          // susceptible being infected
          for (int n=0; n<sus_ages.size(); n++) {
            ppe_multiplier = 1 + (wearing_ppe_susceptible[n] * (ppe_transmission_multiplier - 1));
            lambda[n] += lambda_age[sus_ages[n] - 1] * ppe_multiplier * infection_efficacy[n];
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
