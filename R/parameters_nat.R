#' @title Make vector of variant fold reduction in natural antibody titre values
#' @description Generate a vector for each time step where no VFR effect corresponds
#' to values in the vector of `1`, and the appropriate VFR adjustment otherwise.
#' @param parameters a list of model parameters
#' @param dt time step size
#' @param vfr a positive value giving the fold reduction
#' @param vfr_time_1 time to begin reduction
#' @param vfr_time_2 time to end reduction
#' @export
variant_fold_reduction_vector <- function(parameters, dt, vfr, vfr_time_1, vfr_time_2) {

  stopifnot(is.finite(dt))

  stopifnot(is.finite(vfr))
  stopifnot(is.finite(vfr_time_1))
  stopifnot(is.finite(vfr_time_2))

  stopifnot(vfr > 0)
  stopifnot(vfr_time_1 > 0)
  stopifnot(vfr_time_2 > 0)
  stopifnot(vfr_time_2 > vfr_time_1)
  stopifnot(vfr_time_2 < parameters$time_period)

  vfr_vec <- c(rep(1, times = vfr_time_1 - 1), seq(from = 1, to = vfr, length.out = (vfr_time_2 - vfr_time_1 + 1)), rep(vfr, times = parameters$time_period - vfr_time_2))

  return(vfr_vec)
}


#' @title Attach parameters dealing with immune (NAT) dynamics to parameter list
#' @param parameters a list of model parameters
#' @param vfr a vector from [safir3::variant_fold_reduction_vector]
#' @param mu_ab_infection a vector of values associated with NAT for each infection bout;
#' please read [safir3::attach_event_listeners_natural_immunity] carefully depending on
#' if you are using the additive or overwriting NAT model, because this parameter's
#' interpretation will change. It may also be given as a vector, with each column
#' indexing a day in the simulation (so it is equal to `time_period`).
#' @param std10_infection a standard deviation for log-normal draw of NAT values
#' @export
make_immune_parameters <- function(parameters, vfr, mu_ab_infection = NULL, std10_infection = NULL) {

  if (!is.null(std10_infection)) {
    stopifnot(length(std10_infection) == 1L)
    stopifnot(is.finite(std10_infection))
    stopifnot(std10_infection > 0)
  }

  if (!is.null(mu_ab_infection)) {
    stopifnot(is.finite(mu_ab_infection))
    stopifnot(mu_ab_infection > 0)
    if (inherits(mu_ab_infection, "matrix")) {
      stopifnot(nrow(mu_ab_infection) >= 1)
      stopifnot(ncol(mu_ab_infection) == parameters$time_period)
    } else {
      stopifnot(length(mu_ab_infection) >= 1)
    }
  }

  stopifnot(length(vfr) == parameters$time_period)

  if (is.null(mu_ab_infection)) {
    stopifnot(!is.null(parameters$mu_ab))
    mu_ab_infection <- parameters$mu_ab
  }

  if (is.null(std10_infection)) {
    stopifnot(!is.null(parameters$std10))
    std10_infection <- parameters$std10
  }

  parameters$vfr <- vfr
  parameters$mu_ab_infection <- mu_ab_infection
  parameters$std10_infection <- std10_infection

  return(parameters)
}


#' @title Attach parameters for modeling seperate vaccine and infection derived NAT
#' @param parameters a list of model parameters
#' @param dr_vec_doses a matrix where each column gives the decay rate vector of NATs following a
#' particular dose, and each row is a day
#' @param dr_vec_inf a vector where each element is the decay rate of NATs following infection
#' @param max_ab_inf maximum allowable NAT boost from natural infection
#' @export
make_independent_vaccine_infection_nat_parameters <- function(parameters, dr_vec_doses = NULL, dr_vec_inf = NULL, max_ab_inf = NULL) {
  # check decay for vaccine doses derived NAT
  if (!is.null(dr_vec_doses)) {
    parameters$dr_vec_doses <- dr_vec_doses
  } else {
    stopifnot(!is.null(parameters$dr_vec_doses))
  }
  stopifnot(inherits(parameters$dr_vec_doses, "matrix"))
  stopifnot(nrow(parameters$dr_vec_doses) > 1)
  stopifnot(ncol(parameters$dr_vec_doses) == parameters$max_dose)
  stopifnot(is.finite(parameters$dr_vec_doses))

  # check decay for natural infection derived NAT
  if (!is.null(dr_vec_inf)) {
    parameters$dr_vec_inf <- dr_vec_inf
  } else {
    stopifnot(!is.null(parameters$dr_vec_inf))
  }
  stopifnot(length(parameters$dr_vec_inf) > 1)
  stopifnot(is.finite(parameters$dr_vec_inf))

  # max NAT for infection-derived NAT
  if (!is.null(max_ab_inf)) {
    parameters$max_ab_inf <- max_ab_inf
  } else {
    stopifnot(!is.null(parameters$max_ab_inf))
  }
  stopifnot(is.finite(parameters$max_ab_inf))
  stopifnot(parameters$max_ab_inf > 0)

  return(parameters)
}
