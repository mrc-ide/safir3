# --------------------------------------------------
#   events to handle modeling of natural immunity
#   Sean L. Wu (slwood89@gmail.com)
#   October 2021
# --------------------------------------------------

#' @title Attach event listeners for natural immunity events
#' @param variables a named list of variables
#' @param events a named list of events
#' @param parameters the parameters
#' @param dt size of time step
#' @param additive if `FALSE` the antibody titre is overwritten upon recovery from infection
#' according to `parameters$mu_ab_infection`, if `TRUE` the antibody titre is sampled and added
#' to existing titre. Please choose this option with caution, as it changes the interpretation
#' of the parameter value `parameters$mu_ab_infection`.
#' @export
attach_event_listeners_natural_immunity <- function(variables, events, parameters, dt, additive = FALSE) {

  # checks
  stopifnot(is.logical(additive))

  stopifnot(!is.null(parameters$mu_ab_infection))
  stopifnot(is.finite(parameters$mu_ab_infection))

  if (is.null(parameters$std10_infection)) {
    std10_infection <- parameters$std10
  } else {
    stopifnot(length(parameters$std10_infection) == 1L)
    stopifnot(is.finite(parameters$std10_infection))
    std10_infection <- parameters$std10_infection
  }

  # get current mu_ab_inf
  if (inherits(parameters$mu_ab_infection, "matrix")) {
    get_mu_ab_inf <- function(day) {
      return(parameters$mu_ab_infection[, day])
    }
  } else {
    get_mu_ab_inf <- function(day) {
      return(parameters$mu_ab_infection)
    }
  }

  # effect of infection on NAT
  if (additive) {
    events$recovery$add_listener(
      function(timestep, target) {

        day <- ceiling(timestep * dt)
        mu_ab_inf <- get_mu_ab_inf(day)

        # update inf_num
        inf <- variables$inf_num$get_values(target) + 1L
        variables$inf_num$queue_update(values = inf, index = target)

        # update last time of infection
        variables$inf_time$queue_update(values = timestep, index = target)

        # get NAT values and convert to linear scale
        current_ab_titre <- variables$ab_titre$get_values(index = target)
        current_ab_titre <- exp(current_ab_titre)

        # draw NAT boost on linear scale
        inf[inf > length(mu_ab_inf)] <- length(mu_ab_inf)
        zdose <- 10^rnorm(n = target$size(), mean = log10(mu_ab_inf[inf]),sd = std10_infection)
        new_ab_titre <- current_ab_titre + zdose

        # back to ln scale, and impose max value constraint
        new_ab_titre <- log(new_ab_titre)
        new_ab_titre <- pmin(new_ab_titre, parameters$max_ab)

        # queue NAT update
        variables$ab_titre$queue_update(values = new_ab_titre, index = target)
      }
    )
  } else {
    events$recovery$add_listener(
      function(timestep, target) {

        day <- ceiling(timestep * dt)
        mu_ab_inf <- get_mu_ab_inf(day)

        # update inf_num
        inf <- variables$inf_num$get_values(target) + 1L
        variables$inf_num$queue_update(values = inf, index = target)
        # draw ab titre value
        inf[inf > length(mu_ab_inf)] <- length(mu_ab_inf)
        zdose <- log(10^rnorm(n = target$size(), mean = log10(mu_ab_inf[inf]),sd = std10_infection))
        zdose <- pmin(zdose, parameters$max_ab)
        variables$ab_titre$queue_update(values = zdose, index = target)
        # update last time of infection
        variables$inf_time$queue_update(values = timestep, index = target)
      }
    )
  }

}


# model where infection and vaccine derived NAT stored separately

#' @title Attach event listeners for modeling independent infection-derived NAT
#' @param variables a named list of variables
#' @param events a named list of events
#' @param parameters the parameters
#' @param dt size of time step
#' @export
attach_event_listeners_independent_nat <- function(variables, events, parameters, dt) {

  stopifnot(c("ab_titre_inf", "ab_titre") %in% names(variables))
  stopifnot("max_ab_inf" %in% names(parameters))

  stopifnot(!is.null(parameters$mu_ab_infection))
  stopifnot(is.finite(parameters$mu_ab_infection))

  if (is.null(parameters$std10_infection)) {
    std10_infection <- parameters$std10
  } else {
    stopifnot(length(parameters$std10_infection) == 1L)
    stopifnot(is.finite(parameters$std10_infection))
    std10_infection <- parameters$std10_infection
  }

  # get current mu_ab_inf
  if (inherits(parameters$mu_ab_infection, "matrix")) {
    get_mu_ab_inf <- function(day) {
      return(parameters$mu_ab_infection[, day])
    }
  } else {
    get_mu_ab_inf <- function(day) {
      return(parameters$mu_ab_infection)
    }
  }

  # effect of infection on infection-derived NAT
  events$recovery$add_listener(
    function(timestep, target) {

      day <- ceiling(timestep * dt)
      mu_ab_inf <- get_mu_ab_inf(day)

      # update inf_num
      inf <- variables$inf_num$get_values(target) + 1L
      variables$inf_num$queue_update(values = inf, index = target)

      # update last time of infection
      variables$inf_time$queue_update(values = timestep, index = target)

      # get NAT values and convert to linear scale
      current_nat <- variables$ab_titre_inf$get_values(index = target)
      current_nat <- exp(current_nat)

      # draw NAT boost on linear scale
      inf[inf > length(mu_ab_inf)] <- length(mu_ab_inf)
      nat_boost <- 10^rnorm(n = target$size(), mean = log10(mu_ab_inf[inf]),sd = std10_infection)
      new_nat <- current_nat + nat_boost

      # back to ln scale, and impose max value constraint
      new_nat <- log(new_nat)
      new_nat <- pmin(new_nat, parameters$max_ab_inf)

      # queue NAT update
      variables$ab_titre_inf$queue_update(values = new_nat, index = target)
    }
  )

}

