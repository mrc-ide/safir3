#' @title Create natural immunity variables
#' @param variables a list
#' @param parameters list of model parameters
#' @importFrom individual IntegerVariable
#' @return named list of individual::Variable
#' @export
create_natural_immunity_variables <- function(variables, parameters) {

  n <- parameters$population_size_total

  variables$inf_num <- IntegerVariable$new(initial_values = rep(0L, n))
  variables$inf_time <- IntegerVariable$new(initial_values = rep(-1L, n))

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
