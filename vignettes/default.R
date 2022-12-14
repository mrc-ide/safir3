
library(devtools)
devtools::load_all("C:/xxx/safir3")

# library(safir3)
library(squire)
library(nimue)
library(individual)
library(data.table)
library(ggplot2)
library(dplyr)

start_date <- "2022-11-01"
start_date <- as.Date(start_date)

initial_state <- readRDS("C:/xxx/safir3/data/synthetic_populations/SGP.Rds")
iso3c <- "SGP"
time_period <- 365
dt <- 0.5
initial_infections <- 1000
# scale <- 1 / 1e2

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
  initial_state = initial_state,
  iso3c = iso3c,
  scale = NULL,
  time_period = time_period,
  dt = dt,
  initial_infections = initial_infections
)
timesteps <- parameters$time_period/dt

# vaccine parameters
ab_parameters <- get_vaccine_ab_titre_parameters(vaccine = "Pfizer", max_dose = vaccine_doses,correlated = FALSE)

# combine parameters and verify
parameters <- make_vaccine_parameters(
  safir3_parameters = parameters,
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
hosp_renderer <- Render$new(parameters$time_period)
nat_renderer <- Render$new(parameters$time_period)
dose_renderer <- Render$new(parameters$time_period)
incidence_renderer <- Render$new(timesteps)

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

attach_tracking_listener_incidence(events = events, renderer = incidence_renderer)

# processes
processes <- list(
  natural_immunity_ab_titre_process(parameters = parameters,variables = variables,dt = dt),
  vaccination_process(parameters = parameters,variables = variables,events = events,dt = dt),
  infection_process_vaccine_cpp(parameters = parameters,variables = variables,events = events,dt = dt),
  categorical_count_renderer_process_daily(renderer = renderer,variable = variables$states,categories = variables$states$get_categories(),dt = dt),
  categorical_count_renderer_process_daily(renderer = hosp_renderer,variable = variables$hospitalization_states,categories = variables$hospitalization_states$get_categories(),dt = dt),
  double_count_render_process_daily(renderer = nat_renderer, variable = variables$ab_titre, dt = dt),
  integer_count_render_process_daily(renderer = dose_renderer,variable = variables$dose_num,margin = 0:vaccine_doses,dt = dt)
)

setup_events(parameters = parameters,events = events,variables = variables,dt = dt)

system.time(simulation_loop_safir3(
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

### Health states

saf_dt <- as.data.table(renderer$to_dataframe())
saf_dt[, I_count := I_count]
saf_dt <- melt(saf_dt,id.vars = c("timestep"),variable.name = "name")
saf_dt[, name := gsub("(^)(\\w*)(_count)", "\\2", name)]
setnames(x = saf_dt,old = c("timestep","name","value"),new = c("t","compartment","y"))
saf_dt[["y"]] <- as.integer(saf_dt[["y"]] / parameters$scale)

ggplot(data = saf_dt, aes(t,y,color = compartment)) +
  geom_line() +
  geom_line() +
  facet_wrap(~compartment, scales = "free")

### Hospital states

hosp_saf_dt <- as.data.table(hosp_renderer$to_dataframe())
hosp_saf_dt[, Not_In_Hospital_count := NULL]
hosp_saf_dt[, Pre_Hospital_count := NULL]
hosp_saf_dt <- melt(hosp_saf_dt,id.vars = c("timestep"),variable.name = "name")
hosp_saf_dt[, name := gsub("(^)(\\w*)(_count)", "\\2", name)]
setnames(x = hosp_saf_dt,old = c("timestep","name","value"),new = c("t","compartment","y"))
hosp_saf_dt[["y"]] <- as.integer(hosp_saf_dt[["y"]] / parameters$scale)

ggplot(data = hosp_saf_dt, aes(t,y,color = compartment)) +
  geom_line() +
  geom_line() +
  facet_wrap(~compartment, scales = "free")

### Incidence

df <- incidence_renderer$to_dataframe()
time_steps_per_day = 1 / dt
daily <- data.frame(timestep = rep(0, parameters$time_period), incidence = rep(0, parameters$time_period))
for (t in seq_len(time_steps_per_day)) {
  daily <- daily + df[seq(t, nrow(df), time_steps_per_day), ]
}
daily_incidence <- data.frame(dates = seq.Date(from=start_date, by=1, length.out = parameters$time_period))
daily_incidence$I <- daily$incidence / parameters$scale
daily_incidence[is.na(daily_incidence)] <- 0

library(incidence)
plot(as.incidence(daily_incidence$I, dates = daily_incidence$dates))

# write.csv(daily_incidence,"C:\\tmp\\incidence_safir.csv", row.names = FALSE)
# 
# final_state <- data.frame(age = variables$age$get_values())
# final_state['ab_titre'] <- variables$ab_titre$get_values()
# final_state['infection_number'] <- variables$inf_num$get_values()
# final_state['dose_number'] <- variables$dose_num$get_values()
# final_state['days_since_last_infection'] <- (parameters$time_period - (variables$inf_time$get_values() * dt))
# final_state['days_since_last_dose'] <- (parameters$time_period - (variables$dose_time$get_values() * dt))
# 
# write.csv(final_state,"C:\\tmp\\final_state_safir.csv", row.names = FALSE)

### R estimate

library(EpiEstim)

inf_dist <- parameters$infectiousness_profile
inf_dist[1] <- 0 # EpiEstim requires that si is positive, with value 0 assigned to 1...
si_distr <- inf_dist / sum(inf_dist)
res_non_parametric_si <- estimate_R(daily_incidence, 
                                    method="non_parametric_si",
                                    config = make_config(list(si_distr = si_distr))
)
plot(res_non_parametric_si, "R")