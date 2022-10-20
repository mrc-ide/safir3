
library(devtools)

devtools::load_all("C:/Users/jt511/Documents/GitHub/safir_globalfund")

# library(safir)
library(squire)
library(individual)
library(data.table)
library(ggplot2)
library(dplyr)
require(scales)

iso3c <- "GBR"
scale <- 1 / 1e3
time_period <- 100
dt <- 0.5

# Create parameters
parameters <- get_parameters(
  iso3c = iso3c,
  scale = scale,
  time_period = time_period,
  dt = dt
)
timesteps <- parameters$time_period/dt

# Create variables
variables <- create_variables(parameters = parameters)
variables <- create_infectiousness_variables(variables = variables, parameters = parameters)
variables <- create_hospitalization_variables(variables = variables, parameters = parameters)

# Create events and listeners
events <- create_events(parameters = parameters)
attach_event_listeners(variables = variables,events = events,parameters = parameters, dt = dt)
setup_events(parameters = parameters,events = events,variables = variables,dt = dt)

# Create renderer
renderer <- Render$new(parameters$time_period)

# Create processes
processes <- list(
  infection_process_cpp(parameters = parameters,variables = variables,events = events,dt = dt),
  categorical_count_renderer_process_daily(renderer = renderer,variable = variables$states,categories = variables$states$get_categories(),dt = dt)
  # categorical_count_renderer_process_daily(renderer = renderer,variable = variables$hospitalization_states,categories = variables$hospitalization_states$get_categories(),dt = dt)
)

# Simulation

simulation_loop_safir(
  variables = variables,
  events = events,
  processes = processes,
  timesteps = timesteps,
  variables_dont_update = c("discrete_age"),
  progress = FALSE
)


# Visualization
saf_dt <- as.data.table(renderer$to_dataframe())
saf_dt <- melt(saf_dt,id.vars = c("timestep"),variable.name = "name")
saf_dt[, name := gsub("(^)(\\w*)(_count)", "\\2", name)]
setnames(x = saf_dt,old = c("timestep","name","value"),new = c("t","compartment","y"))
saf_dt[["y"]] <- as.integer(saf_dt[["y"]] / scale)

ggplot(data = saf_dt, aes(t,y,color = compartment)) +
  geom_line() +
  geom_line() +
  facet_wrap(~compartment, scales = "free") + scale_y_continuous(labels = comma)
