# --------------------------------------------------------------------------------
#   infection process for vaccination model (multiple doses, no types)
#   Sean L. Wu (slwood89@gmail.com)
#   July 2021
# --------------------------------------------------------------------------------


#' @title Get vaccine efficacy and Ab titre parameters
#' @param vaccine which vaccine? should be one of: "Pfizer", "AstraZeneca", "Sinovac", "Moderna" (controls the mean titre associated with each dose)
#' @param max_dose maximum number of doses
#' @param correlated are doses correlated?
#' @param hl_s Half life of antibody decay - short
#' @param hl_l Half life of antibody decay - long
#' @param period_s length of the initial decay rate (shortest half-life)
#' @param t_period_l Time point at which to have switched to longest half-life
#' @param ab_50 titre relative to convalescent required to provide 50% protection from infection, on linear scale
#' @param ab_50_severe titre relative to convalescent required to provide 50% protection from severe disease, on linear scale
#' @param std10 Pooled standard deviation of antibody level on log10 data
#' @param k shape parameter of efficacy curve
#' @param nt_transmission_factor used by [safir::vaccine_efficacy_transmission] to compute the effect of antibody titre on onward transmission
#' @param nt_efficacy_transmission a logical variable, specify if antibody titre should affect onward transmission
#' @param max_ab maximum allowable antibody titre draw (on natural log scale)
#' @param mu_ab_list a data.frame
#' @description Get parameters for vaccine efficacy and antibody titre decay rate.
#' @export
get_vaccine_ab_titre_parameters <- function(
  vaccine, max_dose = 2, correlated = FALSE,
  hl_s = 108, hl_l = 3650,
  period_s = 250, t_period_l = 365,
  ab_50 = 0.2,
  ab_50_severe = 0.03,
  std10 = 0.44,
  k = 2.94,
  nt_transmission_factor = 12,
  nt_efficacy_transmission = FALSE,
  max_ab = 8,
  mu_ab_list = data.frame(name = c("Pfizer", "AstraZeneca", "Sinovac", "Moderna"),
                           mu_ab_d1 = c(13/94, 1/59, 28/164, ((185+273)/2)/321),
                           mu_ab_d2 = c(223/94, 32/59, 28/164, 654/158))
) {

  stopifnot(is.logical(nt_efficacy_transmission))

  stopifnot(max_dose %in%  1:(ncol(mu_ab_list) - 1))
  stopifnot(is.logical(correlated))
  stopifnot(vaccine %in% mu_ab_list[, "name"])

  stopifnot(is.finite(nt_transmission_factor))
  stopifnot(nt_transmission_factor > 0)

  stopifnot(is.finite(max_ab))
  stopifnot(max_ab > 0)

  stopifnot(all(is.finite(c(hl_s, hl_l, period_s, t_period_l, ab_50, ab_50_severe, std10, k))))
  stopifnot(all(c(hl_s, hl_l, period_s, t_period_l, ab_50, ab_50_severe, std10, k) > 0))

  time_to_decay <- t_period_l - period_s # time in days to reach longest half-life
  dr_s <- -log(2)/hl_s # Corresponding decay rate in days for half life above
  dr_l <- -log(2)/hl_l

  dr_vec <- c(
    rep(dr_s, period_s),
    seq(dr_s, dr_l, length.out = time_to_decay),
    dr_l
  )

  mu_ab <- as.numeric(mu_ab_list[mu_ab_list$name == vaccine, -1])

  stopifnot(length(mu_ab) >= max_dose)
  stopifnot(all(is.finite(mu_ab)))

  parameters <- list(
    # necessary
    dr_vec = dr_vec,
    mu_ab = mu_ab,
    std10 = std10,
    ab_50 = ab_50,
    ab_50_severe = ab_50_severe,
    k = k,
    nt_transmission_factor = nt_transmission_factor,
    nt_efficacy_transmission = nt_efficacy_transmission,
    max_ab = max_ab,
    correlated = correlated,
    # other
    hl_s = hl_s,
    hl_l = hl_l,
    period_s = period_s,
    t_period_l = t_period_l
  )
  return(parameters)
}


#' @title Combine and verify vaccine parameters
#' @note If modeling a single dose, `dose_period` must be a vector of length 1 and
#' `next_dose_priority_matrix` may be set to `NULL`.
#' @param safir_parameters a list from \code{\link{get_parameters}}
#' @param vaccine_ab_parameters a list from \code{\link{get_vaccine_ab_titre_parameters}}
#' @param vaccine_set a vector giving the number of doses available each day (not each timestep)
#' @param dose_period a vector giving the minimum delay between doses
#' @param strategy_matrix a vaccine strategy matrix from \code{\link[nimue]{strategy_matrix}}
#' or a list of strategy matrices, one for each vaccine dose
#' @param next_dose_priority_matrix a binary matrix giving age groups prioritized for next dose;
#' it should have one fewer row than the number of doses being given, because on the
#' final allocation phase there will be no future dose to prioritize
#' @description Combine parameters for simulation and verify for correctness.
#' @export
make_vaccine_parameters <- function(safir_parameters, vaccine_ab_parameters, vaccine_set, dose_period, strategy_matrix, next_dose_priority_matrix) {

  parameters <- safir_parameters

  vaccine_doses <- length(dose_period)
  stopifnot(vaccine_doses >= 1)
  if (vaccine_doses > 1) {
    stopifnot(all(is.finite(dose_period[2:vaccine_doses])))
    stopifnot(all(dose_period[2:vaccine_doses] >= 0))
  }

  parameters$N_phase <- vaccine_doses
  parameters$dose_period <- dose_period

  stopifnot(length(vaccine_set) == parameters$time_period)
  stopifnot(all(is.finite(vaccine_set)))
  stopifnot(all(vaccine_set >= 0))

  parameters$vaccine_set <- floor(vaccine_set)

  storage.mode(next_dose_priority_matrix) <- "integer"
  stopifnot(nrow(next_dose_priority_matrix) == vaccine_doses - 1)

  # one strategy matrix for all dose phases (input matrix)
  if (inherits(strategy_matrix, 'matrix')) {

    N_age <- ncol(strategy_matrix)
    strategy_matrix_array <- replicate(n = vaccine_doses, expr = strategy_matrix, simplify = FALSE)

    # one strategy matrix for each dose phase (input list)
  } else if(inherits(strategy_matrix, 'list')) {

    stopifnot(length(strategy_matrix) == vaccine_doses)

    N_age <- ncol(strategy_matrix[[1]])
    strategy_matrix_array <- strategy_matrix

    # bad input
  } else {
    stop("invalid object passed for strategy_matrix")
  }

  stopifnot(ncol(next_dose_priority_matrix) == N_age)
  stopifnot(parameters$number_of_age_groups == N_age)

  parameters$vaccine_coverage_mat <- strategy_matrix_array
  parameters$next_dose_priority <- next_dose_priority_matrix

  # attach necessary parameters from vaccine pars
  parameters$dr_vec <- vaccine_ab_parameters$dr_vec
  parameters$mu_ab <- vaccine_ab_parameters$mu_ab
  parameters$std10 <- vaccine_ab_parameters$std10
  parameters$ab_50 <- vaccine_ab_parameters$ab_50
  parameters$ab_50_severe <- vaccine_ab_parameters$ab_50_severe
  parameters$k <- vaccine_ab_parameters$k
  parameters$nt_transmission_factor <- vaccine_ab_parameters$nt_transmission_factor
  parameters$nt_efficacy_transmission <- vaccine_ab_parameters$nt_efficacy_transmission
  parameters$max_ab <- vaccine_ab_parameters$max_ab
  parameters$correlated <- vaccine_ab_parameters$correlated

  return(parameters)
}
