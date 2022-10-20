
library(devtools)

devtools::load_all("C:/Users/jthompson/Documents/GitHub/safir_globalfund")

# library(safir)
library(squire)
library(safir)
library(nimue)
library(individual)
library(data.table)
library(ggplot2)
library(dplyr)

iso3c <- "GBR"
scale <- 1 / 1e3
time_period <- 365
dt <- 0.5

# vaccine dosing
vaccine_doses <- 2
dose_period <- c(NaN, 28)
vaccine_set <- c(0, seq(from = 1e3, to = 1e4, length.out = time_period-1))
vaccine_set <- floor(vaccine_set)

# vaccine strategy
vaccine_coverage_mat <- strategy_matrix(strategy = "Elderly",max_coverage = 0.0)
next_dose_priority <- matrix(data = 0, nrow = vaccine_doses - 1,ncol = ncol(vaccine_coverage_mat))
next_dose_priority[1, 15:17] <- 1 # prioritize 3 oldest age groups for next dose

# base parameters
parameters <- get_parameters(
  iso3c = iso3c,
  scale = scale,
  time_period = time_period,
  dt = dt
)
timesteps <- parameters$time_period/dt

# vaccine parameters
ab_parameters <- get_vaccine_ab_titre_parameters(vaccine = "Pfizer", max_dose = vaccine_doses,correlated = FALSE)

# combine parameters and verify
parameters <- make_vaccine_parameters(
  safir_parameters = parameters,
  vaccine_ab_parameters = ab_parameters,
  vaccine_set = vaccine_set,
  dose_period = dose_period,
  strategy_matrix = vaccine_coverage_mat,
  next_dose_priority_matrix = next_dose_priority
)

# ab boost for each infection
parameters$mu_ab_infection <- ab_parameters$mu_ab

# create variables
variables <- create_variables(parameters = parameters)
variables <- create_infectiousness_variables(variables = variables, parameters = parameters)
variables <- create_hospitalization_variables(variables = variables, parameters = parameters)
variables <- create_natural_immunity_variables(variables = variables, parameters = parameters)
variables <- create_vaccine_variables(variables = variables, parameters = parameters)

# create events
events <- create_events(parameters = parameters)
events <- create_events_vaccination(events = events, parameters = parameters)
attach_event_listeners(variables = variables, events = events, parameters = parameters, dt = dt)
attach_event_listeners_vaccination(variables = variables, events = events, parameters = parameters, dt = dt)
attach_event_listeners_natural_immunity(variables = variables, events = events, parameters = parameters, dt = dt)

# make renderers
renderer <- Render$new(parameters$time_period)
nat_renderer <- Render$new(parameters$time_period)
dose_renderer <- Render$new(parameters$time_period)

double_count_render_process_daily <- function(renderer, variable, dt) {
  stopifnot(inherits(variable, "DoubleVariable"))
  stopifnot(inherits(renderer, "Render"))
  function(t) {
    if ((t * dt) %% 1 == 0) {
      day <- as.integer(t * dt)
      nat <- exp(variable$get_values())
      quantiles <- quantile(x = nat, probs = c(0.025, 0.5, 0.975))
      renderer$render(name = "q025", value = quantiles[[1]], timestep = day)
      renderer$render(name = "q5", value = quantiles[[2]], timestep = day)
      renderer$render(name = "q975", value = quantiles[[3]], timestep = day)
      renderer$render(name = "mean", value = mean(x = nat), timestep = day)
    }
  }
}

# processes
processes <- list(
  natural_immunity_ab_titre_process(parameters = parameters,variables = variables,dt = dt),
  vaccination_process(parameters = parameters,variables = variables,events = events,dt = dt),
  infection_process_vaccine_cpp(parameters = parameters,variables = variables,events = events,dt = dt),
  categorical_count_renderer_process_daily(renderer = renderer,variable = variables$states,categories = variables$states$get_categories(),dt = dt),
  double_count_render_process_daily(renderer = nat_renderer, variable = variables$ab_titre, dt = dt),
  integer_count_render_process_daily(renderer = dose_renderer,variable = variables$dose_num,margin = 0:vaccine_doses,dt = dt)
)

setup_events(parameters = parameters,events = events,variables = variables,dt = dt)

system.time(simulation_loop_safir(
  variables = variables,
  events = events,
  processes = processes,
  timesteps = timesteps,
  variables_dont_update = c("discrete_age", "phase"),
  progress = FALSE
))

## Plot Results

### Antibody titre

ab_titre_dt <- as.data.table(nat_renderer$to_dataframe())
setnames(ab_titre_dt, "timestep", "Day")

ggplot(data = ab_titre_dt) +
  geom_line(aes(x=Day,y=mean)) +
  geom_ribbon(aes(x=Day,ymin=q025,ymax=q975),alpha=0.5) +
  theme_bw()

### Proportion of population with each dose

dose_out <- dose_renderer$to_dataframe()
colnames(dose_out)[2:(vaccine_doses+2)] <- as.character(0:vaccine_doses)
dose_out <- melt(as.data.table(dose_out),id.vars="timestep")
setnames(dose_out, "variable", "dose")

ggplot(data = dose_out) +
  geom_line(aes(x=timestep,y=value,color=dose)) +
  theme_bw()

### Infection states

saf_dt <- as.data.table(renderer$to_dataframe())
saf_dt[, I_count := I_count]
saf_dt <- melt(saf_dt,id.vars = c("timestep"),variable.name = "name")
saf_dt[, name := gsub("(^)(\\w*)(_count)", "\\2", name)]
setnames(x = saf_dt,old = c("timestep","name","value"),new = c("t","compartment","y"))
saf_dt[["y"]] <- as.integer(saf_dt[["y"]] / scale)

ggplot(data = saf_dt, aes(t,y,color = compartment)) +
  geom_line() +
  geom_line() +
  facet_wrap(~compartment, scales = "free")
