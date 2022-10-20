# --------------------------------------------------
#   functions to output data for analysis
#   Sean L. Wu (slwood89@gmail.com)
#   June 2021
# --------------------------------------------------

#' @title Create listener to track incidence
#'
#' @param renderer a \code{\link[individual]{Render}} object
#' @export
create_incidence_tracking_listener <- function(renderer) {
  stopifnot(inherits(renderer, "Render"))
  function(timestep, target) {
    renderer$render(name = "incidence", value = target$size(), timestep = timestep)
  }
}

#' @title Attach listener to track incidence
#' @description To use incidence tracking, run the lines in examples
#' anytime after all events have been made and initialized and before
#' simulation begins.
#' @param events a list of events
#' @param renderer a \code{\link[individual]{Render}} object
#' @examples
#' \dontrun{
#' incidence_renderer <- individual::Render$new(timesteps)
#' attach_tracking_listener_incidence(events = events, renderer = incidence_renderer)
#' }
#' @export
attach_tracking_listener_incidence <- function(events, renderer) {
  events$infection$add_listener(
    create_incidence_tracking_listener(renderer = renderer)
  )
}

#' @title Create listener to track age-stratified incidence
#'
#' @param renderer a [individual::Render] object
#' @param age an [individual::IntegerVariable] object
#' @param parameters model parameters ([safir::get_parameters])
#' @export
create_age_incidence_tracking_listener <- function(renderer, age, parameters) {
  stopifnot(inherits(parameters, "list"))
  stopifnot(inherits(renderer, "Render"))
  stopifnot(inherits(age, "IntegerVariable"))
  N_age <- parameters$number_of_age_groups
  function(timestep, target) {
    for (a in seq_len(N_age)) {
      age_bset <- age$get_index_of(a)
      age_bset$and(target)
      renderer$render(name = paste0("incidence_age_", a), value = age_bset$size(), timestep = timestep)
    }
  }
}

#' @title Attach listener to track age-stratified incidence
#' @description To use age-stratified incidence tracking, run the lines in examples
#' anytime after all events have been made and initialized but before
#' simulation begins.
#' @param events a list of events
#' @param renderer a [individual::Render] object
#' @param age an [individual::IntegerVariable] object
#' @param parameters model parameters ([safir::get_parameters])
#' @examples
#' \dontrun{
#' incidence_renderer <- individual::Render$new(timesteps)
#' attach_tracking_listener_age_incidence(
#'     events = events, renderer = age_incidence_renderer,
#'     age = variables$discrete_age, parameters = parameters
#'  )
#' }
#' @export
attach_tracking_listener_age_incidence <- function(events, renderer, age, parameters) {
  events$infection$add_listener(
    create_age_incidence_tracking_listener(renderer = renderer, age = age, parameters = parameters)
  )
}


#' @title Render integer variables
#' @description Renders the number of individuals in each integer count. Should only be used
#' for integer variables with bounded ranges
#' @param renderer a \code{\link[individual]{Render}} object
#' @param variable a \code{\link[individual]{IntegerVariable}} object
#' @param margin the set of values which may be taken for the integer variable
#' @export
integer_count_render_process <- function(renderer, variable, margin) {
  stopifnot(inherits(variable, "IntegerVariable"))
  stopifnot(inherits(renderer, "Render"))
  function(t) {
    for (c in margin) {
      renderer$render(paste0(c, '_count'), variable$get_size_of(c), t)
    }
  }
}

#' @title Render integer variables every day
#' @description Renders the number of individuals in each integer count. Should only be used
#' for integer variables with bounded ranges. This only renders output on timesteps
#' that correspond to a day. This saves memory and potentially time for
#' simulations with small timesteps.
#' @param renderer a \code{\link[individual]{Render}} object
#' @param variable a \code{\link[individual]{IntegerVariable}} object
#' @param margin the set of values which may be taken for the integer variable
#' @param dt size of time step
#' @export
integer_count_render_process_daily <- function(renderer, variable, margin, dt) {
  stopifnot(inherits(variable, "IntegerVariable"))
  stopifnot(inherits(renderer, "Render"))
  function(t) {
    if ((t * dt) %% 1 == 0) {
      for (c in margin) {
        renderer$render(paste0(c, '_count'), variable$get_size_of(c), as.integer(t * dt))
      }
    }
  }
}

#' @title Render doses by age every day
#' @description Renders the number of individuals with each dose in each age bin
#' (cross tabulated doses and ages).
#' This only renders output on timesteps that correspond to a day.
#' See examples for how to quickly summarize the output.
#' @param renderer a \code{\link[individual]{Render}} object
#' @param age a \code{\link[individual]{IntegerVariable}} object
#' @param dose a \code{\link[individual]{IntegerVariable}} object
#' @param parameters model parameters
#' @param dt size of time step
#' @examples
#' \dontrun{
#' # if the renderer object is called dose_age_renderer
#' tmp <- as.data.table(dose_age_renderer$to_dataframe())
#' tmp <- melt(tmp,id.vars="timestep")
#' tmp1 <- tmp[, tstrsplit(variable, "_", keep = c(2, 4))]
#' tmp[ , variable := NULL]
#' tmp <- cbind(tmp, tmp1)
#' setnames(x = tmp,old = c("V1", "V2"),new = c("dose","age"))
#' dose_dt <- tmp[, .(value = sum(value)), by = .(dose, timestep)]
#' }
#' @export
dose_age_render_process_daily <- function(renderer, age, dose, parameters, dt) {
  stopifnot(inherits(age, "IntegerVariable"))
  stopifnot(inherits(dose, "IntegerVariable"))
  stopifnot(inherits(renderer, "Render"))
  num_dose <- parameters$N_phase
  num_age <- parameters$number_of_age_groups
  function(t) {
    if ((t * dt) %% 1 == 0) {
      day <- as.integer(t * dt)
      dose_age_tab <- cross_tab_doses_age(doses = dose$.variable, age = age$.variable, num_doses = num_dose, num_ages = num_age)
      for (d in seq_len(nrow(dose_age_tab))) {
        for (a in seq_len(ncol(dose_age_tab))) {
          renderer$render(paste0("dose_", d-1, "_age_", a, "_count"), dose_age_tab[d, a], day)
        }
      }
    }
  }
}

#' @title Render categories every day
#' @description Renders the number of individuals in each category, but only
#' on timesteps that correspond to a day. This saves memory and potentially time
#' for simulations with small timesteps.
#' @param renderer a \code{\link[individual]{Render}} object
#' @param variable a \code{\link[individual]{CategoricalVariable}} object
#' @param categories a character vector of categories to render
#' @param dt size of time step
#' @return a function which can be passed as a process to \code{\link{simulation_loop}}
#' @export
categorical_count_renderer_process_daily <- function(renderer, variable, categories, dt) {
  stopifnot(inherits(variable, "CategoricalVariable"))
  stopifnot(inherits(renderer, "Render"))
  function(t) {
    if ((t * dt) %% 1 == 0) {
      for (c in categories) {
        renderer$render(paste0(c, '_count'), variable$get_size_of(c), as.integer(t * dt))
      }
    }
  }
}


#' @title Render categories by age every day
#' @description This only renders output on timesteps that correspond to a day.
#' See examples for how to quickly summarize the output.
#' @param renderer a [individual::Render] object
#' @param age a [individual::IntegerVariable] object
#' @param compartments a [individual::CategoricalVariable] object
#' @param parameters model parameters
#' @param dt size of time step
#' @examples
#' \dontrun{
#' # if the renderer object is called comp_render
#' tmp <- as.data.table(comp_render$to_dataframe())
#' tmp <- melt(tmp,id.vars="timestep")
#' tmp1 <- tmp[, tstrsplit(variable, "_", keep = c(2, 4))]
#' tmp[ , variable := NULL]
#' tmp <- cbind(tmp, tmp1)
#' setnames(x = tmp,old = c("V1", "V2"),new = c("compartment","age"))
#' compartment_dt <- tmp[, .(value = sum(value)), by = .(compartment, timestep)]
#' }
#' @export
compartments_age_render_process_daily <- function(renderer, age, compartments, parameters, dt) {
  stopifnot(inherits(age, "IntegerVariable"))
  stopifnot(inherits(compartments, "CategoricalVariable"))
  stopifnot(inherits(renderer, "Render"))
  num_age <- parameters$number_of_age_groups
  compartment_names <- compartments$get_categories()
  function(t) {
    if ((t * dt) %% 1 == 0) {
      day <- as.integer(t * dt)
      comp_age_tab <- cross_tab_compartments_age(compartments = compartments$.variable, age = age$.variable, compartment_names = compartment_names, num_ages = num_age)
      for (c in seq_len(ncol(comp_age_tab))) {
        for (a in seq_len(nrow(comp_age_tab))) {
          renderer$render(paste0("compartment_", compartment_names[c], "_age_", a, "_count"), comp_age_tab[a, c], day)
        }
      }
    }
  }
}
