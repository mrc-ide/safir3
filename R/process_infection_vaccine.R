# -------------------------------------------------------
#   infection process for individual infectiousness model
# -------------------------------------------------------

#' @title C++ infection process (modified squire transmission model)
#'
#' @description Simulates the infection process for the squire transmission model, modified
#'              to include variable individual infectiousness and antiviral treatment
#' Calls \code{\link{infection_process_vaccine_cpp_internal}} to return an external pointer object.
#'
#' @param parameters Model parameters
#' @param variables a list of model variables, the output of [create_variables]
#' @param events a list of [individual::TargetedEvent], the output of [create_events]
#' @param dt the time step
#' @export
infection_process_vaccine_cpp <- function(parameters, variables, events, dt) {

  stopifnot(all(c("states","discrete_age") %in% names(variables)))
  stopifnot("infection" %in% names(events))

  return(
    infection_process_vaccine_cpp_internal(
      parameters = parameters,
      variables = variables,
      states = variables$states$.variable,
      discrete_age = variables$discrete_age$.variable,
      infectiousness_start_time = variables$infectiousness_start_time$.variable,
      wearing_ppe = variables$wearing_ppe$.variable,
      infection = events$infection$.event,
      dt = dt
    )
  )
}
