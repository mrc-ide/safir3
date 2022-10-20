# --------------------------------------------------
#   vaccination process for the general vaccine model
#   but with no types
#   Sean L. Wu (slwood89@gmail.com)
#   June 2021
# --------------------------------------------------


#' @title Vaccination distribution process (multiple doses, no types)
#' @description This function distributes vaccine doses each day for the vaccine
#' model with multiple doses but no types. It distributes vaccines according to the
#' following strategy.
#'
#'   * start at \code{phase = 1}, meaning we begin by distributing the first dose
#'   * use \code{\link{get_vaccination_priority_stage}} to check what stage of the
#'     prioritization matrix \code{vaccine_coverage_mat} from \code{\link[nimue]{strategy_matrix}}
#'     we should be using to decide who to vaccinate.
#'   * we can advance one row of the prioritization matrix when coverage for this dose (phase) >=
#'     prioritization matrix and coverage for next dose (phase) >= prioritization matrix groups prioritized
#'     for the next dose by the matrix \code{next_dose_priority}.
#'
#' @param parameters a named list
#' @param variables a named list
#' @param events a named list
#' @param dt size of time step
#' @export
vaccination_process <- function(parameters, variables, events, dt) {

  stopifnot(all(c("dose_num","dose_time","phase","discrete_age") %in% names(variables)))
  stopifnot(all(c("N_phase", "dose_period", "vaccine_coverage_mat", "next_dose_priority", "vaccine_set") %in% names(parameters)))
  stopifnot(is.finite(parameters$N_phase))
  stopifnot(nrow(parameters$next_dose_priority) == parameters$N_phase - 1)
  stopifnot(ncol(parameters$next_dose_priority) == ncol(parameters$vaccine_coverage_mat[[1]]))

  return(

    function(timestep) {

      day <- timestep * dt

      # only distribute once a day
      if (day %% 1 == 0) {

        # get phase we are on
        phase <- variables$phase$value

        # how many doses we have today
        doses_left <- parameters$vaccine_set[day]

        # if vaccination done or no doses available, return early
        if (phase > parameters$N_phase || doses_left <= 0) {
          return(invisible(NULL))
        }

        # calculate what step we are on
        step <- get_vaccination_priority_stage(variables = variables, events = events, phase = phase, parameters = parameters)

        # advance to next phase: recalculate (step = -1 is a signal to advance to the next dosing phase)
        while (step == -1) {
          variables$phase$value <- variables$phase$value + 1L
          phase <- variables$phase$value
          # if advancing phase goes over the total number, vaccination is done and we return
          if (phase > parameters$N_phase) {
            return(invisible(NULL))
          }
          step <- get_vaccination_priority_stage(variables = variables, events = events, phase = phase, parameters = parameters)
        }

        stopifnot(phase <= parameters$N_phase)

        if (parameters$N_phase > 1 & phase < parameters$N_phase) {

          priority_ages <- parameters$next_dose_priority[phase, ]
          if (sum(priority_ages) > 0) {

            N_prioritisation_steps_priority <- nrow(parameters$vaccine_coverage_mat[[phase + 1L]])

            # step through priority matrix for *next* dose
            step_priority <- 0L

            while (doses_left > 0 & step_priority < N_prioritisation_steps_priority) {

              step_priority <- step_priority + 1L
              p_step <- parameters$vaccine_coverage_mat[[phase + 1L]][step_priority, ]

              eligible <- target_pop(
                dose = phase + 1L, variables = variables, events = events, parameters = parameters,
                timestep = timestep, dt = dt, strategy_matrix_step = p_step, next_dose_priority = priority_ages
              )
              doses_left <- assign_doses(
                doses_left = doses_left, events = events, dose = phase + 1L,
                eligible = eligible, parameters = parameters
              )

            }

          }

        }

        # row of the vaccine coverage matrix
        vaccine_coverage_mat <- parameters$vaccine_coverage_mat[[phase]]
        N_prioritisation_steps <- nrow(vaccine_coverage_mat)
        p_step <- vaccine_coverage_mat[step, ]

        # get eligible persons for this dose phase and strategy step
        eligible <- target_pop(
          dose = phase, variables = variables, events = events, parameters = parameters,
          timestep = timestep, dt = dt, strategy_matrix_step = p_step
        )
        # assign doses to eligible persons based on available supply
        doses_left <- assign_doses(
          doses_left = doses_left, events = events, dose = phase,
          eligible = eligible, parameters = parameters
        )

        # if remaining doses, give out for this dose phase according to the prioritization matrix
        while (doses_left > 0 & step < N_prioritisation_steps) {

          step <- step + 1
          p_step <- vaccine_coverage_mat[step, ]

          eligible <- target_pop(
            dose = phase, variables = variables, events = events, parameters = parameters,
            timestep = timestep, dt = dt, strategy_matrix_step = p_step
          )
          doses_left <- assign_doses(
            doses_left = doses_left, events = events, dose = phase,
            eligible = eligible, parameters = parameters
          )

        }

      } # end daily vaccine distribution

    } # end function

  ) # end return

}
