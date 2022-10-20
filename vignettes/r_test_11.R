
library(devtools)

devtools::load_all("C:/Users/jt511/Documents/GitHub/safir_globalfund")

# library(safir)
library(squire)
library(individual)
library(data.table)
library(ggplot2)
library(dplyr)
library(GA)
require(scales)
library(zoo)

time_period_fortnights <- 10
iso3c <- "GBR"
scale <- 1 / 1e3
time_period <- 14 * time_period_fortnights
dt <- 0.5

# Set target
target <- rep(2000, time_period)

# Create parameters
parameters <- get_parameters(
  iso3c = iso3c,
  scale = scale,
  time_period = time_period,
  dt = dt
)
timesteps <- parameters$time_period/dt

run <- function(basic_beta_set) {

  # Process beta_set input
  beta_set = rep(basic_beta_set, each=14)
  beta_set <- rollmean(beta_set, k = 7)
  first_val <- beta_set[1]
  last_val <- beta_set[length(beta_set)]
  beta_set <- c(first_val, first_val, first_val, beta_set, last_val, last_val, last_val)
  parameters$beta_set <- beta_set

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
    categorical_count_renderer_process_daily(renderer = renderer,variable = variables$states,categories = c("D"),dt = dt)
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

  # Process deaths output
  saf_dt <- as.data.table(renderer$to_dataframe())
  output <- c(0, diff(as.integer(saf_dt[["D_count"]] / scale)))
  first_val <- output[1]
  last_val <- output[length(output)]
  output <- c(first_val, first_val, first_val, rollmean(output, k = 7), last_val, last_val, last_val)

  # Calculate fitness
  1 / (sum((output - target) ^ 2))

}

time_period_fortnights <- 10
GA <- ga(type = "real-valued",
         fitness = function(x) run(x),
         lower = rep(0.0, time_period_fortnights), upper = rep(1.0, time_period_fortnights),
         popSize = 50, maxiter = 100, run = 20)

summary(GA)
plot(GA)

best_sol <- GA@solution
