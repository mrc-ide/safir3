#' @title Create natural immunity variables
#' @param variables a list
#' @param parameters list of model parameters
#' @importFrom individual IntegerVariable
#' @return named list of individual::Variable
#' @export
create_natural_immunity_variables <- function(variables, parameters) {

  if(is.null(parameters$initial_state)){
    n <- parameters$population_size_total
    variables$inf_num <- IntegerVariable$new(initial_values = rep(0L, n))
    variables$inf_time <- IntegerVariable$new(initial_values = rep(-1L, n))
  } else {
    infection_number <- parameters$initial_state[['infection_number']]
    days_since_last_infection <- parameters$initial_state[['days_since_last_infection']]
    # days_since_last_infection[days_since_last_infection > 900] <- -1
    ticks_since_last_infection <- as.integer(days_since_last_infection / parameters$dt)
    variables$inf_num <- IntegerVariable$new(initial_values = infection_number)
    variables$inf_time <- IntegerVariable$new(initial_values = -ticks_since_last_infection)
  }

  return(variables)
}


#' @title Create independent infection-derived NAT variables
#' @description This creates an independent set of variables to track infection-derived NAT independently
#' from vaccine-derived NAT.
#' @param variables a list
#' @param parameters list of model parameters
#' @importFrom individual DoubleVariable
#' @return named list of individual::Variable
#' @export
create_independent_nat_variables <- function(variables, parameters) {

  n <- parameters$population_size_total

  variables <- create_natural_immunity_variables(variables = variables, parameters = parameters)

  # ab dynamics
  variables$ab_titre_inf <- DoubleVariable$new(initial_values = rep(-Inf, n))

  return(variables)
}
