# --------------------------------------------------
#   functions to work with scheduling/allocating vaccines
#   for the general vaccine model but with no types
#   Sean L. Wu (slwood89@gmail.com)
#   June 2021
# --------------------------------------------------


#' @title Get people who either have this dose or are scheduled for it
#' @description
#' Return a list of bitsets, one for each age group, containing.
#' those persons who either have had this dose or are scheduled for it.
#' The complement of each of these sets will be those
#' people who haven't gotten this dose and aren't scheduled for it.
#' @param variables a list of model variables
#' @param events a list of model events
#' @param dose what dose to check coverage for
#' @param parameters a list of parameters
#' @export
get_current_coverage <- function(variables, events, dose, parameters) {

  bsets <- vector("list", parameters$number_of_age_groups)

  vaccinated_bset <- variables$dose_num$get_index_of(set = dose:parameters$N_phase) # have gotten this dose
  vaccinated_bset$or(events$scheduled_dose[[dose]]$get_scheduled()) # have gotten this dose OR are scheduled for it

  for (a in 1:parameters$number_of_age_groups) {

    age_bset <- variables$discrete_age$get_index_of(a)
    vaccinated_age_bset <- vaccinated_bset$copy()
    vaccinated_age_bset$and(age_bset) # in the right age

    bsets[[a]] <- vaccinated_age_bset

  }

  return(bsets)
}


#' @title Get persons who are eligible for this dose
#' @description Get people who are eligible for this dose, from the output of \code{\link{get_current_coverage}}.
#' Because \code{coverage} is people (by age) who either have this dose or are scheduled for it,
#' if we are checking eligibility for dose 1 we just want the complement of those people. For doses > 1 we need to
#' first take the complement and then check
#' to see who of those persons got their previous dose past the waiting threshold.
#' This function returns a list of bitsets.
#' @param timestep the current time step
#' @param dt time step size
#' @param coverage output of \code{\link{get_current_coverage}}
#' @param variables a list of model variables
#' @param dose what dose to get eligibility for
#' @param parameters a list of parameters
#' @importFrom individual Bitset
#' @export
get_current_eligible_from_coverage <- function(timestep, dt, coverage, variables, dose, parameters) {

  # eligible persons are those not covered and in this age group
  bsets <- lapply(X = 1:parameters$number_of_age_groups,FUN = function(a){
    not_cov <- coverage[[a]]$copy()
    not_cov$not(inplace = TRUE)
    not_cov$and(variables$discrete_age$get_index_of(a))
  })

  # for doses > 1 we need to check if they are past threshold
  if (dose > 1) {

    threshold <- timestep - as.integer(parameters$dose_period[dose]/dt)

    # nobody eligible if timestep is too soon for threshold to have elapsed
    if (threshold < 0) {

      for (a in 1:parameters$number_of_age_groups) {
        bsets[[a]] <- Bitset$new(size = bsets[[a]]$max_size)
      }

    # threshold past 0
    } else {

      previous_dose_in_threshold <- variables$dose_num$get_index_of(set = dose - 1)
      previous_dose_in_threshold$and(variables$dose_time$get_index_of(a = 0, b = threshold))
      for (a in 1:parameters$number_of_age_groups) {
        bsets[[a]]$and(previous_dose_in_threshold)
      }

    }

  }

  # for the 1st dose, we don't need to check previous threshold, we can just return

  return(bsets)
}


#' @title Get prioritisation step for a specific dosing phase
#' @description Return which row of the \code{\link[nimue]{strategy_matrix}} the vaccination
#' program should be targeting for coverage. To complete a step, all \code{phase}
#' dose coverage should be > prioritisation matrix target and all \code{phase + 1}
#' dose coverage for prioritized groups should be > prioritisation matrix target.
#' If these two conditions are fulfilled for the entire \code{phase}, the function returns
#' -1 to indicate that vaccination dosing \code{phase} should be advanced.
#' @param variables a list of model variables
#' @param events a list of model events
#' @param phase what dosing phase are we on
#' @param parameters a list of parameters
#' @export
get_vaccination_priority_stage <- function(variables, events, phase, parameters) {

  # age groups
  age_size <- parameters$population_size_by_age_group

  # strategy matrix for this phase
  vaccine_coverage_mat <- parameters$vaccine_coverage_mat[[phase]]
  N_prioritisation_steps <- nrow(vaccine_coverage_mat)

  stopifnot(is.finite(parameters$N_phase))

  # calculate coverage for this dose
  coverage_this_dose <- get_current_coverage(variables = variables, events = events, dose = phase, parameters = parameters)

  pr_this_dose <- get_size_bset(bsets = coverage_this_dose)
  pr_this_dose <- pr_this_dose / age_size

  # go through prioritization steps
  for (p in 1:N_prioritisation_steps) {

    ages_2_check <- which(vaccine_coverage_mat[p, ] > 0)
    this_dose_not_cover <- any(pr_this_dose[ages_2_check] < vaccine_coverage_mat[p, ages_2_check])

    if (this_dose_not_cover) {
      return(p)
    }
  }

  # if we did not return by now it means this step is complete, return -1
  return(-1)

}


#' @title Target persons in each age group to vaccinate (multi-dose, no types)
#' @description Given the current dosing phase (\code{dose}), the current step of the strategy matrix (\code{strategy_matrix_step})
#' and possible the age groups prioritized for the next dose \code{next_dose_priority}, figure out who should get this \code{dose}.
#' This will return a list of bitsets giving people who are eligible to get this dose, taking into account that they
#' have not yet had it, nor are scheduled for it, and (if dose > 1) are past the threshold for their previous dose.
#' @param dose what dose to target (the phase we are vaccinating)
#' @param variables a list of model variables
#' @param events a list of model events
#' @param parameters a list of parameters
#' @param timestep the current time step
#' @param dt the size of the time step
#' @param strategy_matrix_step the current step of the strategy matrix
#' @param next_dose_priority (optional) a row of the next dose priority matrix
#' @export
target_pop <- function(dose, variables, events, parameters, timestep, dt, strategy_matrix_step, next_dose_priority = NULL) {

  # age groups
  age_size <- parameters$population_size_by_age_group

  # get people who have gotten this dose or are scheduled for it
  coverage_this_dose <- get_current_coverage(
    variables = variables, events = events, dose = dose, parameters = parameters
  )

  # current coverage proportion
  current_coverage <- get_size_bset(bsets = coverage_this_dose)
  current_coverage <- current_coverage / age_size

  # who is eligible for this dose?
  eligible_this_dose <- get_current_eligible_from_coverage(
    timestep = timestep, dt = dt, coverage = coverage_this_dose, variables = variables, dose = dose, parameters = parameters
  )

  # coverage targets are the strategy matrix step; potentially multiplied by next_dose_priority
  coverage_targets <- strategy_matrix_step

  if (!is.null(next_dose_priority)) {
    coverage_targets <- coverage_targets * next_dose_priority
  }

  # remaining population left to cover with current dose number (dose) to reach target coverage in prioritization step
  # this does not take into account who is eligible based on threshold
  n_to_cover <- ceiling(pmax(0, (coverage_targets - current_coverage)) * age_size)

  # final number to target for vaccine is minimum of eligible[a] and n_to_cover[a]
  eligible_num <- get_size_bset(bsets = eligible_this_dose)
  n_to_cover <- pmin(n_to_cover, eligible_num)

  # if nobody to target, set empty bitset, otherwise use choose to randomly select
  for (a in 1:parameters$number_of_age_groups) {
    if (n_to_cover[a] < 1) {
      eligible_this_dose[[a]] <- Bitset$new(size = eligible_this_dose[[a]]$max_size)
    } else {
      if (n_to_cover[a] < eligible_this_dose[[a]]$size()) {
        eligible_this_dose[[a]]$choose(k = n_to_cover[a])
      }
      stopifnot(n_to_cover[a] <= eligible_this_dose[[a]]$size()) # assertion
    }
  }

  return(eligible_this_dose)
}


#' @title Assign doses to eligible persons based on available supply (multi-dose, no types)
#' @description Given a list of bitsets indicating who is eligible for this dose from \code{\link{target_pop}},
#' and the number of doses available today \code{doses_left}, schedule individuals for doses.
#' If fewer doses are available than eligible persons, the doses are apportioned
#' proportional to age group size.
#' @param doses_left the number of vaccine doses available
#' @param events a list of model events
#' @param dose what dose to assign for
#' @param eligible output of \code{\link{target_pop}}
#' @param parameters a list of model parameters
#' @importFrom stats rmultinom
#' @export
assign_doses <- function(doses_left, events, dose, eligible, parameters) {

  n_to_assign <- get_size_bset(bsets = eligible)

  # check if can quit early
  if (sum(n_to_assign) < 1) {
    return(doses_left)
  }

  # no dose scarcity: everyone can get one
  if (sum(n_to_assign) <= doses_left) {

    for (a in 1:parameters$number_of_age_groups) {

      if (n_to_assign[a] < 1) {
        next()
      } else {
        events$scheduled_dose[[dose]]$schedule(target = eligible[[a]], delay = 0)
        doses_left <- doses_left - n_to_assign[a]
      }

    }

  # dose scarcity: allocate proportional to group size
  } else {

    # assign limited doses to age groups
    group_weights <- n_to_assign / sum(n_to_assign)
    assigned <- floor(doses_left * group_weights)

    if(sum(assigned) != doses_left){
      extra_to_assign <- n_to_assign - assigned
      extra_weights <- extra_to_assign / sum(extra_to_assign)
      assigned <- assigned + as.vector(rmultinom(n = 1,size = doses_left - sum(assigned),prob = extra_weights))
    }

    stopifnot(sum(assigned) <= doses_left)

    for (a in 1:parameters$number_of_age_groups) {

      if (assigned[a] < 1) {
        next()
      } else {
        eligible[[a]]$choose(k = assigned[a])
        stopifnot(eligible[[a]]$size() == assigned[a])
        events$scheduled_dose[[dose]]$schedule(target = eligible[[a]], delay = 0)
        doses_left <- doses_left - assigned[a]
      }

    }

  } # end dose scarcity

  return(doses_left)
}
