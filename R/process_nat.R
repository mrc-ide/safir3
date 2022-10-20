# update NAT for vaccine-derived NAT only

#' @title Process that updates the antibody (Ab) titre each time step
#' @description The values in `ab_titre` are calculated on the log scale.
#' This process will not calculate decay correctly for `dt > 1` so that is disallowed.
#' @param parameters a list of model parameters
#' @param variables a list of model variables
#' @param dt time step size
#' @export
vaccine_ab_titre_process <- function(parameters, variables, dt) {

  stopifnot(dt <= 1)

  return(
    function(timestep) {

      # only those with at least 1 dose will need to have titre calculated
      vaccinated <- variables$dose_num$get_index_of(set = 0)
      vaccinated$not(inplace = TRUE)

      if (vaccinated$size() > 0) {

        # for each person we need to know the time since their last dose
        time_since_last_dose <- get_time_since_last_dose(timestep = timestep, dt = dt, vaccinated = vaccinated, dose_time = variables$dose_time)

        # ceiling to go to next integer day and do not exceed the end of decay rate vector
        time_since_last_dose <- ceiling(time_since_last_dose)
        time_since_last_dose[time_since_last_dose > length(parameters$dr_vec)] <- length(parameters$dr_vec)

        # current Ab titre
        current_ab_titre <- variables$ab_titre$get_values(index = vaccinated)

        # new Ab titre
        new_ab_titre <- current_ab_titre + (parameters$dr_vec[time_since_last_dose] * dt)

        # schedule an update
        variables$ab_titre$queue_update(values = new_ab_titre, index = vaccinated)

      }

    }
  )

}

#' @title Calculate the time elapsed in days since each person's last dose
#' @param timestep current time step
#' @param dt size of time step
#' @param vaccinated bitset of vaccinated persons
#' @param dose_time a list of \code{\link[individual]{IntegerVariable}} objects
#' @export
get_time_since_last_dose <- function(timestep, dt, vaccinated, dose_time) {

  # output just for the vaccinated persons
  times <- rep(NaN, vaccinated$size())
  times <- timestep - dose_time$get_values(index = vaccinated)
  times <- times * dt

  return(times)
}


# update NAT for combined vaccine and infection derived NAT

#' @title Process that updates the antibody (Ab) titre each time step for model vaccine model with natural immunity
#' @description The values in `ab_titre` are calculated on the log scale.
#' This process will not calculate decay correctly for `dt > 1` so that is disallowed.
#' @param parameters a list of model parameters
#' @param variables a list of model variables
#' @param dt time step size
#' @export
natural_immunity_ab_titre_process <- function(parameters, variables, dt) {

  return(
    function(timestep) {

      # persons who have been vaccinated or infected will have antibody titre
      vaccinated <- variables$dose_num$get_index_of(set = 0)$not(inplace = TRUE)
      infected <- variables$inf_num$get_index_of(set = 0)$not(inplace = TRUE)

      vaccinated_or_infected <- vaccinated$or(infected)

      if (vaccinated_or_infected$size() > 0) {

        # for each person we need to know the time since their last dose
        time_since_last_dose_or_infection <- get_time_since_last_dose_or_infection(
          timestep = timestep, dt = dt, vaccinated_or_infected = vaccinated_or_infected,
          dose_time = variables$dose_time, inf_time = variables$inf_time
        )

        # ceiling to go to next integer day and do not exceed the end of decay rate vector
        time_since_last_dose_or_infection <- ceiling(time_since_last_dose_or_infection)
        time_since_last_dose_or_infection[time_since_last_dose_or_infection > length(parameters$dr_vec)] <- length(parameters$dr_vec)

        # current Ab titre
        current_ab_titre <- variables$ab_titre$get_values(index = vaccinated_or_infected)

        # new Ab titre
        new_ab_titre <- current_ab_titre + (parameters$dr_vec[time_since_last_dose_or_infection] * dt)

        # schedule an update
        variables$ab_titre$queue_update(values = new_ab_titre, index = vaccinated_or_infected)

      }

    }
  )

}

#' @title Calculate the time elapsed in days since each person's last dose or infection
#' @param timestep current time step
#' @param dt size of time step
#' @param vaccinated_or_infected [individual::Bitset] of vaccinated persons
#' @param dose_time an [individual::IntegerVariable] object
#' @param inf_time an [individual::IntegerVariable] object
#' @export
get_time_since_last_dose_or_infection <- function(timestep, dt, vaccinated_or_infected, dose_time, inf_time) {

  last_dose_tt <- dose_time$get_values(vaccinated_or_infected)
  last_inf_tt <- inf_time$get_values(vaccinated_or_infected)

  # the max of last dose and last inf time will give the time of last ab boost
  times <- timestep - pmax(last_dose_tt, last_inf_tt)
  times <- times * dt

  return(times)
}


# update NAT for independent vaccine and infection derived NAT

#' @title Update NAT for seperate vaccine and infection-derived NAT
#' @description Process that updates NAT titre for model where vaccine and infection derived
#' NAT are stored separately.
#' The NAT values are stored on the natural log scale.
#' This process will not calculate decay correctly for `dt > 1` so that is disallowed.
#' @param parameters a list of model parameters
#' @param variables a list of model variables
#' @param dt time step size
#' @export
independent_ab_titre_process <- function(parameters, variables, dt) {

  stopifnot(c("ab_titre_inf", "ab_titre", "dose_num", "inf_num", "dose_time", "inf_time") %in% names(variables))
  stopifnot("dr_vec_doses" %in% names(parameters))
  stopifnot(inherits(parameters$dr_vec_doses, "matrix"))
  stopifnot("dr_vec_inf" %in% names(parameters))

  stopifnot(inherits(variables$ab_titre_inf, "DoubleVariable"))
  stopifnot(inherits(variables$ab_titre, "DoubleVariable"))
  stopifnot(inherits(variables$dose_num, "IntegerVariable"))
  stopifnot(inherits(variables$inf_num, "IntegerVariable"))
  stopifnot(inherits(variables$dose_time, "IntegerVariable"))
  stopifnot(inherits(variables$inf_time, "IntegerVariable"))

  return(
    function(timestep) {

      # persons who have been vaccinated or infected will have NAT
      vaccinated <- variables$dose_num$get_index_of(set = 0)$not(inplace = TRUE)
      infected <- variables$inf_num$get_index_of(set = 0)$not(inplace = TRUE)

      # update NATs for people who are only vaccinated
      if (vaccinated$size() > 0) {

        # for each person we need to know the time since their last dose
        time_since_last_dose <- get_time_since_last_dose(timestep = timestep, dt = dt, vaccinated = vaccinated, dose_time = variables$dose_time)

        # last dose
        last_dose_num <- variables$dose_num$get_values(index = vaccinated)

        # ceiling to go to next integer day and do not exceed the end of decay rate vector
        time_since_last_dose <- ceiling(time_since_last_dose)
        time_since_last_dose[time_since_last_dose > nrow(parameters$dr_vec_doses)] <- nrow(parameters$dr_vec_doses)

        # current Ab titre
        current_nat <- variables$ab_titre$get_values(index = vaccinated)

        # new Ab titre
        dr_vec_doses_index <- do.call(what = cbind, args = list(time_since_last_dose, last_dose_num))
        new_nat <- current_nat + (parameters$dr_vec_doses[dr_vec_doses_index] * dt)

        # schedule an update
        variables$ab_titre$queue_update(values = new_nat, index = vaccinated)

      }

      # update NATs for people who have only been infected
      if (infected$size() > 0) {

        time_since_last_infection <- get_time_since_last_infection(timestep = timestep, dt = dt, infected = infected, inf_time = variables$inf_time)

        time_since_last_infection <- ceiling(time_since_last_infection)
        time_since_last_infection[time_since_last_infection > length(parameters$dr_vec_inf)] <- length(parameters$dr_vec_inf)

        current_nat <- variables$ab_titre_inf$get_values(index = infected)

        new_nat <- current_nat + (parameters$dr_vec_inf[time_since_last_infection] * dt)

        # schedule an update
        variables$ab_titre_inf$queue_update(values = new_nat, index = infected)

      }

    } # end return function
  )

}


# independent_ab_titre_process <- function(parameters, variables, dt) {
#
#   stopifnot(c("ab_titre_inf", "ab_titre", "dose_num", "inf_num", "dose_time", "inf_time") %in% names(variables))
#   stopifnot("dr_vec_doses" %in% names(parameters))
#   stopifnot(inherits(parameters$dr_vec_doses, "matrix"))
#   stopifnot("dr_vec_inf" %in% names(parameters))
#
#   stopifnot(inherits(variables$ab_titre_inf, "DoubleVariable"))
#   stopifnot(inherits(variables$ab_titre, "DoubleVariable"))
#   stopifnot(inherits(variables$dose_num, "IntegerVariable"))
#   stopifnot(inherits(variables$inf_num, "IntegerVariable"))
#   stopifnot(inherits(variables$dose_time, "IntegerVariable"))
#   stopifnot(inherits(variables$inf_time, "IntegerVariable"))
#
#   return(
#     function(timestep) {
#
#       # persons who have been vaccinated or infected will have NAT
#       vaccinated <- variables$dose_num$get_index_of(set = 0)$not(inplace = TRUE)
#       infected <- variables$inf_num$get_index_of(set = 0)$not(inplace = TRUE)
#
#       # bitsets for vaxx/inf, vaxx only, inf only
#       vaxx_and_inf <- vaccinated$copy()$and(infected)
#       vaccinated$set_difference(vaxx_and_inf)
#       infected$set_difference(vaxx_and_inf)
#
#       # update NATs for people who are only vaccinated
#       if (vaccinated$size() > 0) {
#
#         # for each person we need to know the time since their last dose
#         time_since_last_dose <- get_time_since_last_dose(timestep = timestep, dt = dt, vaccinated = vaccinated, dose_time = variables$dose_time)
#
#         # last dose
#         last_dose_num <- variables$dose_num$get_values(index = vaccinated)
#
#         # ceiling to go to next integer day and do not exceed the end of decay rate vector
#         time_since_last_dose <- ceiling(time_since_last_dose)
#         time_since_last_dose[time_since_last_dose > nrow(parameters$dr_vec_doses)] <- nrow(parameters$dr_vec_doses)
#
#         # current Ab titre
#         current_nat <- variables$ab_titre$get_values(index = vaccinated)
#
#         # new Ab titre
#         dr_vec_doses_index <- do.call(what = cbind, args = list(time_since_last_dose, last_dose_num))
#         new_nat <- current_nat + (parameters$dr_vec_doses[dr_vec_doses_index] * dt)
#
#         # schedule an update
#         variables$ab_titre$queue_update(values = new_nat, index = vaccinated)
#
#       }
#
#       # update NATs for people who have only been infected
#       if (infected$size() > 0) {
#
#         time_since_last_infection <- get_time_since_last_infection(timestep = timestep, dt = dt, infected = infected, inf_time = variables$inf_time)
#
#         time_since_last_infection <- ceiling(time_since_last_infection)
#         time_since_last_infection[time_since_last_infection > length(parameters$dr_vec_inf)] <- length(parameters$dr_vec_inf)
#
#         current_nat <- variables$ab_titre_inf$get_values(index = infected)
#
#         new_nat <- current_nat + (parameters$dr_vec_inf[time_since_last_infection] * dt)
#
#         # schedule an update
#         variables$ab_titre_inf$queue_update(values = new_nat, index = infected)
#
#       }
#
#       # update NATs for people who have been vaccinated and infected
#       if (vaxx_and_inf$size() > 0) {
#
#         # update infection-derived NATs
#         time_since_last_infection <- get_time_since_last_infection(timestep = timestep, dt = dt, infected = vaxx_and_inf, inf_time = variables$inf_time)
#
#         # make sure no times exceeed the end of the respective decay vectors
#         time_since_last_infection <- ceiling(time_since_last_infection)
#         time_since_last_infection[time_since_last_infection > length(parameters$dr_vec_inf)] <- length(parameters$dr_vec_inf)
#
#         current_inf_nat <- variables$ab_titre_inf$get_values(index = vaxx_and_inf)
#
#         new_inf_nat <- current_inf_nat + (parameters$dr_vec_inf[time_since_last_infection] * dt)
#
#         # schedule an update
#         variables$ab_titre_inf$queue_update(values = new_inf_nat, index = vaxx_and_inf)
#
#         # update vaccine-derived NATs
#         time_since_last_dose <- get_time_since_last_dose(timestep = timestep, dt = dt, vaccinated = vaxx_and_inf, dose_time = variables$dose_time)
#
#         # make sure no times exceeed the end of the respective decay vectors
#         time_since_last_dose <- ceiling(time_since_last_dose)
#         time_since_last_dose[time_since_last_dose > nrow(parameters$dr_vec_doses)] <- nrow(parameters$dr_vec_doses)
#
#         # last dose
#         last_dose_num <- variables$dose_num$get_values(index = vaxx_and_inf)
#
#         # current Ab titre
#         current_vaxx_nat <- variables$ab_titre$get_values(index = vaxx_and_inf)
#
#         # new Ab titre
#         dr_vec_doses_index <- do.call(what = cbind, args = list(time_since_last_dose, last_dose_num))
#         new_vaxx_nat <- current_vaxx_nat + (parameters$dr_vec_doses[dr_vec_doses_index] * dt)
#
#         # schedule an update
#         variables$ab_titre$queue_update(values = new_vaxx_nat, index = vaxx_and_inf)
#
#       }
#
#     } # end return function
#   )
#
# }

#' @title Calculate the time elapsed in days since each person's last infection
#' @param timestep current time step
#' @param dt size of time step
#' @param infected [individual::Bitset] of infected persons
#' @param inf_time an [individual::IntegerVariable] object
#' @export
get_time_since_last_infection <- function(timestep, dt, infected, inf_time) {

  last_inf_tt <- inf_time$get_values(infected)

  # the max of last dose and last inf time will give the time of last ab boost
  times <- timestep - last_inf_tt
  times <- times * dt

  return(times)
}

