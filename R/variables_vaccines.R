# --------------------------------------------------
#   Create and query vaccination variables
#   Sean L. Wu (slwood89@gmail.com)
#   June 2021
# --------------------------------------------------


# --------------------------------------------------
#   nimue version
# --------------------------------------------------

#' @title Create vaccination variables (nimue vaccine model)
#' @description Create all individual variables for humans
#'
#' @param variables a list from [create_variables]
#' @param pop population list
#' @importFrom individual Bitset
#' @importFrom individual IntegerVariable
#' @return named list of individual::Variable
#' @export
create_vaccine_variables_nimue <- function(variables, pop) {

  n <- sum(pop$n)

  variables$vaccine_states <- IntegerVariable$new(initial_values = rep(1,n))
  variables$eligible <- Bitset$new(size = n)
  variables$vaccinated <- Bitset$new(size = n)
  variables$empty <- Bitset$new(size = n)

  return(variables)
}


#' @title Get proportion of an age group that is vaccinated (nimue vaccine model)
#' @description Get proportion of an age group that has received a particular vaccine dose
#' by this timestep.
#' @param variables a list
#' @param age an age group (can be a vector of multiple age groups)
#'
#' @export
get_proportion_vaccinated_nimue <- function(variables, age) {
  age_bset <- variables$discrete_age$get_index_of(age)
  N <- age_bset$size()
  vaccinated_bset <- variables$vaccinated$copy()
  vaccinated_bset$and(age_bset)
  return( vaccinated_bset$size() / N )
}


# --------------------------------------------------
#   multiple doses, no types
# --------------------------------------------------

#' @title Create vaccination variables (multi-dose, no types)
#' @param variables a list
#' @param parameters list of model parameters
#' @importFrom individual IntegerVariable
#' @importFrom individual DoubleVariable
#' @return named list of individual::Variable
#' @export
create_vaccine_variables <- function(variables, parameters) {

  n <- parameters$population_size_total
  max_dose <- parameters$N_phase
  correlated <- parameters$correlated

  if(is.null(parameters$initial_state)){

    variables$dose_num <- IntegerVariable$new(initial_values = rep(0,n))
    variables$dose_time <- IntegerVariable$new(initial_values = rep(-1, n))
    variables$phase <- new.env(hash = FALSE)
    variables$phase$value <- 1L

    # ab dynamics
    variables$ab_titre <- DoubleVariable$new(initial_values = rep(-Inf, n))

    if (correlated) {
      variables$zdose <- DoubleVariable$new(initial_values = rep(-Inf, n))
    }

  } else {

    dose_number <- parameters$initial_state[['dose_number']]
    days_since_last_dose <- parameters$initial_state[['days_since_last_dose']]
    ab_titre <- parameters$initial_state[['ab_titre']]
    ticks_since_last_dose <- as.integer(days_since_last_dose / parameters$dt)
    variables$dose_num <- IntegerVariable$new(initial_values = dose_number)
    variables$dose_time <- IntegerVariable$new(initial_values = -ticks_since_last_dose)
    variables$phase <- new.env(hash = FALSE)
    variables$phase$value <- 1L

    # ab dynamics
    variables$ab_titre <- DoubleVariable$new(initial_values = ab_titre)

    if (correlated) {
      variables$zdose <- DoubleVariable$new(initial_values = rep(-Inf, n))
    }

  }

  return(variables)
}


#' @title Initialize vaccination variables (multi-dose, no types)
#' @param variables a list from \code{\link{create_vaccine_variables}}
#' @param dose_time_init list of dose times (must all be integer vectors of the same length, -1 is value to indicate that individual has not
#' yet received this dose)
#' @param dose_num_init vector of initial doses (must be self-consistent with \code{dose_time_init})
#' @export
initialize_vaccine_variables <- function(variables, dose_time_init, dose_num_init) {
  stopifnot(inherits(dose_time_init, "list"))
  stopifnot(length(dose_time_init) > 0)
  stopifnot(length(dose_time_init[[1]]) == length(dose_num_init))
  # if you haven't gotten the first dose, you are at dose 0
  stopifnot(all(which(dose_time_init[[1]] < 0) == which(dose_num_init == 0)))
  for (d in seq_along(dose_time_init)) {
    dose_times <- dose_time_init[[d]]
    variables$dose_time$queue_update(values = dose_times[which(dose_times >= 0)], index = which(dose_times >= 0))
    variables$dose_time$.update()
  }
  variables$dose_num$queue_update(value = dose_num_init)
  variables$dose_num$.update()
}


#' @title Update vaccine variables
#' @description This should be called from the simulation loop. It does not
#' update disease state.
#' @param variables a list
#' @export
update_vaccine_variables <- function(variables) {

  variables$dose_time$.update()
  variables$dose_num$.update()
  variables$ab_titre$.update()

  if (!is.null(variables$zdose)) {
    variables$zdose$.update()
  }

}
