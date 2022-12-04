# --------------------------------------------------
#   model parameters
# --------------------------------------------------


#' @title Get parameters
#'
#' @param iso3c Country iso3c
#' @param scale Scale factor
#' @param time_period Length of simulation (in days)
#' @param dt Size of time step
#' @param initial_infections Number of intial infections
#' @param ... Other parameters
#' @return model parameters
#' @export
get_parameters <- function(initial_state = NULL,
                           iso3c = NULL,
                           scale = 1.0,
                           time_period = 365,
                           dt = 0.5,
                           initial_infections = 1000,
                           ...) {

    stopifnot(is.finite(dt) & dt > 0)
    stopifnot(!is.null(iso3c))

    pars <- list()

    # iso3c and country name
    pars$iso3c <- iso3c
    pars$country <- squire::population[squire::population$iso3c == iso3c, "country"][[1]]

    if(is.null(initial_state)){

      # Scale factor
      pars$scale <- scale

      # Population data (from Squire)
      assert_string(iso3c)
      if(!iso3c %in% unique(squire::population$iso3c)){stop("iso3c not found")}
      pars$pop <- squire::population[squire::population$iso3c == iso3c, ] %>%
                  dplyr::arrange(.data$age_group)
      pars$pop$n <- as.integer(pars$pop$n * scale)

      # Population size
      pars$population_size_by_age_group <- pars$pop$n
      pars$population_size_total <- sum(pars$pop$n)
      pars$number_of_age_groups <- length(pars$population_size_by_age_group)

      # Max age
      pars$max_age <- 100

    } else {

      # Initial state
      pars$initial_state <- initial_state

      # Population size
      age <- pars$initial_state[['age']]
      pars$population_size_total <- length(age)
      pars$number_of_age_groups <- 17
      population_size_by_age_group <- rep(0, pars$number_of_age_groups)
      for(a in age){
        discrete_age <- min(a %/% 5, 16) + 1
        population_size_by_age_group[discrete_age] <- population_size_by_age_group[discrete_age] + 1
      }
      pars$population_size_by_age_group <- population_size_by_age_group

      # Max age
      pars$max_age <- max(age)

      # Population data (from Squire)
      assert_string(iso3c)
      if(!iso3c %in% unique(squire::population$iso3c)){stop("iso3c not found")}
      pars$pop <- squire::population[squire::population$iso3c == iso3c, ] %>%
                  dplyr::arrange(.data$age_group)

      # Scale factor
      pars$scale <- pars$population_size_total / sum(pars$pop$n)
    }

    # Time period (in days)
    pars$time_period <- time_period

    # Time step (in days)
    pars$dt <- dt

    # Contact matrix (from Squire)
    contact_matrix <- squire::contact_matrices[[pars$pop$matrix[1]]]
    contact_matrix <- process_contact_matrix_scaled_age(contact_matrix,
                                                        pars$population_size_by_age_group)
    contact_matrix <- div_pop(contact_matrix, pars$population_size_by_age_group)
    pars$contact_matrix <- contact_matrix

    # External force of infection
    pars$lambda_external <- rep(0, time_period)

    # Input sequence of betas (issue: mixing matrices not consistent, so this depends on country...)
    pars$beta_set <- rep(0.5, time_period)

    # Transition probabilities by age group (each set should sum to 1)
    pars$prob_IAsymp         <- rep(0.4, 17)
    pars$prob_IMild          <- rep(0.4, 17)
    pars$prob_IModerate      <- rep(0.1, 17)
    pars$prob_ISevere        <- rep(0.085, 17)
    pars$prob_ICritical      <- rep(0.015, 17) # Redundant
    pars$prob_IAsymp_to_S    <- rep(1.0, 17)
    pars$prob_IAsymp_to_D    <- rep(0.0, 17) # Redundant
    pars$prob_IMild_to_S     <- rep(1.0, 17)
    pars$prob_IMild_to_D     <- rep(0.0, 17) # Redundant
    pars$prob_IModerate_to_S <- rep(1.0, 17)
    pars$prob_IModerate_to_D <- rep(0.0, 17) # Redundant
    pars$prob_ISevere_to_S   <- rep(0.95, 17)
    pars$prob_ISevere_to_D   <- rep(0.05, 17) # Redundant
    pars$prob_ICritical_to_S <- rep(0.5, 17)
    pars$prob_ICritical_to_D <- rep(0.5, 17) # Redundant

    # Number of initial infections
    pars$num_initial_infections <- max(as.integer(initial_infections * pars$scale), 1)

    # Initial number in each health state
    pars$S_0 <- pars$population_size_total - pars$num_initial_infections
    pars$I_0 <- pars$num_initial_infections
    pars$D_0 <- 0

    # Infectiousness (durations need not coincide with the length of the infectiousness profile)
    pars$infectiousness_profile <- c(0.000395985, 0.001768802, 0.125221560, 0.222678871,
                                     0.210222459, 0.176880156, 0.125221560, 0.070417258,
                                     0.039598534, 0.017688016, 0.007041726, 0.002226789,
                                     0.000994670)
    pars$dur_IAsymp    <- 12 / pars$dt
    pars$dur_IMild     <- 12 / pars$dt
    pars$dur_IModerate <- 12 / pars$dt
    pars$dur_ISevere   <- 12 / pars$dt
    pars$dur_ICritical <- 12 / pars$dt

    # Hospitalization
    pars$prob_hosp_IAsymp    <- 0.0
    pars$prob_hosp_IMild     <- 0.0
    pars$prob_hosp_IModerate <- 0.0
    pars$prob_hosp_ISevere   <- 0.5
    pars$prob_hosp_ICritical <- 0.9
    pars$rr_hosp_IAsymp      <- 0.5 # risk reduction to be calculated here using decision tree
    pars$rr_hosp_IMild       <- 0.5 # risk reduction to be calculated here using decision tree
    pars$rr_hosp_IModerate   <- 0.5 # risk reduction to be calculated here using decision tree
    pars$rr_hosp_ISevere     <- 0.5 # risk reduction to be calculated here using decision tree
    pars$rr_hosp_ICritical   <- 0.5 # risk reduction to be calculated here using decision tree
    pars$enter_hosp_delay    <- 5 / pars$dt # time between infection and entering hospital
    pars$dur_hospital        <- 12 / pars$dt # time between entering and leaving hospital

    # Hospital capacity
    pars$beds_total <- get_healthcare_capacity(country = pars$country,
                                               pop_size = pars$population_size_total)

    # PPE
    pars$ppe_transmission_multiplier <- 0.1

    return(pars)

}


#' Divide a contact matrix by population
#'
#' @param contact_matrix Matrix
#' @param population Population vector
#'
#' @return Matrix
div_pop <- function(contact_matrix, population){
  t(t(contact_matrix) / population)
}


#' Process a contact matrix
#'
#' @param contact_matrix A contact matrix
#' @param population Vector of population by age
#'
#' @return Processed matrix
process_contact_matrix_scaled_age <- function(contact_matrix, population) {

  # Convert unbalanced matrix of per-capita rates to total number of contacts
  # between different age groups and balance by taking the mean of i->j and j->i
  contact_matrix <- rbind(contact_matrix, contact_matrix[16,])
  contact_matrix <- cbind(contact_matrix, contact_matrix[,16]*population[17] / sum(population[16:17]))
  contact_matrix[,16] <- contact_matrix[,16]*population[16] / sum(population[16:17])
  MIJ <- t(vapply(seq(population),function(x){
    contact_matrix[x,] * population[x]
  }, FUN.VALUE = numeric(length(population))))
  adjust_mat <- (MIJ + t(MIJ))/2 # symmetric and balanced

  # Convert to new per-capita rates by dividing by Population
  # Resulting Matrix Is Asymmetric But Balanced
  # Asymmetric in that c_ij != c_ji BUT Total Number of Contacts i->j and j->i
  # Is Balanced (so when we divide by pop at end, will be balanced)
  processed_matrix <- t(vapply(seq(population), function(x) {
    adjust_mat[x, ] / population[x]
  }, FUN.VALUE = numeric(length(population))))

  # Adjusting to create input for model i.e. per capita rates divided by
  # population to give the number of contacts made on each individual
  return(processed_matrix)
}


#' Set hospital bed capacity
#'
#' @param country Country name
#' @param pop_size Total number of people in the population
#'
#' @return Number of hospital beds for this country
get_healthcare_capacity <-  function(country, pop_size){

  if(country %in% unique(squire::country_specific_healthcare_capacity$country)) {
    beds <- squire::country_specific_healthcare_capacity[match(country,
                                            squire::country_specific_healthcare_capacity$country), ]
  } else {
    income_group <- squire::income_group$income_group[match(country, squire::income_group$country)]
    beds <- squire::income_strata_healthcare_capacity[
                           squire::income_strata_healthcare_capacity$income_group == income_group, ]
  }
  beds_total <- as.integer(beds$hosp_beds * pop_size / 1000)

  return(beds_total)
}