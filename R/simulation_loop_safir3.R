# --------------------------------------------------
#   Simulation loop for nimue style vaccinations
#   Sean L. Wu (slwood89@gmail.com)
#   May 2021
# --------------------------------------------------

#' @title Simulation loop for safir3 models
#' @description This is a replacement for [individual::simulation_loop]
#' for models where not all values in `variables` need to be updated, and
#' also adds an optional progress bar
#' @param variables a list of Variables
#' @param events a list of Events
#' @param processes a list of processes to execute on each timestep
#' @param timesteps the number of timesteps to simulate
#' @param variables_dont_update character vector of variables that won't be updated
#' @param progress dispaly a progress bar?
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @export
simulation_loop_safir3 <- function(
  variables = list(),
  events = list(),
  processes = list(),
  timesteps,
  variables_dont_update,
  progress = FALSE
) {

  stopifnot(is.finite(timesteps))
  stopifnot(timesteps > 0)

  to_update <- which(!(names(variables) %in% variables_dont_update))
  stopifnot(length(to_update) > 0)

  if (progress) {
    pb <- txtProgressBar(min = 1, max = timesteps, initial = 1, style = 3)
  }

  to_process <- which(names(events) != "scheduled_dose")
  stopifnot(length(to_process) > 0)

  # sim loop
  for (t in seq_len(timesteps)) {

    # execute processes
    for (process in processes) {
      execute_any_process(process, t)
    }

    # fire events
    for (event in events[to_process]) {
      event$.process()
    }
    if (!is.null(events$scheduled_dose)) {
      for (event in events$scheduled_dose) {
        event$.process()
      }
    }

    # update state
    for (variable in variables[to_update]) {
      variable$.update()
    }

    # event clocks tick
    for (event in events[to_process]) {
      event$.tick()
    }
    if (!is.null(events$scheduled_dose)) {
      for (event in events$scheduled_dose) {
        event$.tick()
      }
    }

    # progress bar
    if (progress) {
      setTxtProgressBar(pb = pb,value = t)
    }
  }
  # end loop, close out pb
  if (progress) {
    close(con = pb)
  }
}
