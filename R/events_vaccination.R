# --------------------------------------------------
#   create event functions full vaccine model
#   Sean L. Wu (slwood89@gmail.com)
#   June 2021
# --------------------------------------------------

# write a function

#' @title Append vaccination events
#'
#' @param events a named list of individual::Event
#' @param parameters model parameters
#' @importFrom individual TargetedEvent
#' @export
create_events_vaccination <- function(events, parameters) {

  # pop size
  N <- parameters$population_size_total

  # scheduled future doses
  events$scheduled_dose <- replicate(n = parameters$N_phase,expr = TargetedEvent$new(N),simplify = FALSE)

  return(events)
}


#' @title Schedule some individuals for a vaccination dose
#' @description This is called from the event listeners for each dose, and also
#' aids in better testing of the simulation model
#' @param timestep current time step
#' @param variables a list
#' @param target a \code{\link[individual]{Bitset}}
#' @param dose which dose
#' @param parameters model parameters
#' @importFrom stats rnorm
#' @export
schedule_dose_vaccine <- function(timestep, variables, target, dose, parameters) {

  variables$dose_num$queue_update(value = dose,index = target)
  variables$dose_time$queue_update(values = timestep, index = target)

  if (inherits(target,"Bitset")) {
    n <- target$size()
  } else {
    n <- length(target)
  }

  nat <- variables$ab_titre$get_values(index = target) # get the stored NAT for that individual
  nat <- exp(nat) # transform back to linear scale

  if (parameters$correlated) {

    if (dose > 1) {

      # correlated doses > 1; use ratio of mean titre
      zdose_prev <- variables$zdose$get_values(index = target) # get the drawn titre value from previous dose (this is on natural log scale)

      zdose_prev_linear <- exp(zdose_prev)

      zdose_linear <- (zdose_prev_linear * (parameters$mu_ab[dose] / parameters$mu_ab[dose-1L])) # get value of new dose based on previous dose as they are correlated, linear scale

      # check it doesnt exceed our maximum allowed titre
      zdose_linear <- pmin(zdose_linear, exp(parameters$max_ab))

      # store the zdose value on nat log scale
      variables$zdose$queue_update(values = log(zdose_linear), index = target)

      # update nat based on difference between dose i and dose i - 1
      nat <- nat + zdose_linear

      # transform back to natural log scale
      nat <- log(nat)

      # check new nat doesnt exceed max allowed titre
      nat <- pmin(nat, parameters$max_ab)

    } else {

      # initial dose titre on linear scale
      zdose <- 10^rnorm(n = n, mean = log10(parameters$mu_ab[1L]), sd = parameters$std10)

      # check zdose doesnt exceed max value
      zdose <- pmin(zdose, exp(parameters$max_ab))

      # store the zdose value on nat log scale
      variables$zdose$queue_update(values = log(zdose), index = target)

      # update NAT
      nat <- nat + zdose

      # transform back to natural log scale
      nat <- log(nat)

      # check new nat doesnt exceed max allowed titre
      nat <- pmin(nat, parameters$max_ab)

    }

  } else {

    zdose_linear <- 10^rnorm(n = n, mean = log10(parameters$mu_ab[dose]), sd = parameters$std10) # get value of new dose

    # check it doesnt exceed our maximum allowed titre
    zdose_linear <- pmin(zdose_linear, exp(parameters$max_ab))

    # update nat based on difference between dose i and dose i - 1
    nat <- nat + zdose_linear

    # transform back to natural log scale
    nat <- log(nat)

    # check new nat doesnt exceed max allowed titre
    nat <- pmin(nat, parameters$max_ab)

  }

  variables$ab_titre$queue_update(values = nat, index = target)

}


#' @title Create listener for vaccination dose (multi-dose, no types)
#' @description Updates state when a vaccine dose is given. It does not schedule future events.
#' @param variables a named list of variables
#' @param parameters the parameters
#' @param dose integer dose
create_vaccination_dose_listener <- function(variables, parameters, dose) {
  stopifnot( dose > 0 )
  function (timestep, target) {
    schedule_dose_vaccine(timestep = timestep,variables = variables,target = target,dose = dose, parameters = parameters)
  }
}


#' @title Attach event listeners for vaccination events (multi-dose, no types)
#' @param variables a named list of variables
#' @param events a named list of events
#' @param parameters the parameters
#' @param dt size of time step
#' @export
attach_event_listeners_vaccination <- function(variables, events, parameters, dt) {

  # add to each dose
  for (d in seq_along(events$scheduled_dose)) {
    events$scheduled_dose[[d]]$add_listener(
      create_vaccination_dose_listener(variables = variables,parameters = parameters,dose = d)
    )
  }

}
