# --------------------------------------------------
#   event functions
# --------------------------------------------------


#' @title Create events
#' @description Creates a list of all events in the model
#' @param parameters model parameters
#' @export
create_events <- function(parameters) {

    N <- parameters$population_size_total

    list(
        infection = TargetedEvent$new(N),
        recovery = TargetedEvent$new(N),
        death = TargetedEvent$new(N),
        enter_hospital_waiting = TargetedEvent$new(N),
        enter_hospital = TargetedEvent$new(N),
        leave_hospital = TargetedEvent$new(N)
    )

}


#' @title Determines the outcome of infections
#' @description Schedules hosptialization, recovery and death events
#' @param prob_recover Baseline probability of recovery (no hospitalization)
#' @param prob_hosp Probability of hospitalization
#' @param rr_hosp Risk reduction of death associated to hospitalization
#' @param target A bitset of newly infected individuals
#' @param availability Available hospital beds
#' @param R_delay A function representing the distribution of time to recovery
#' @param D_delay A function representing the distribution of time to death
#' @param availability Available hospital beds
#' @export
determine_outcome <- function(prob_recover, prob_hosp, rr_hosp, target, availability, S_delay,
                              D_delay, events, parameters, calculate_nat, day) {

    hosp <- target$copy()

    # NAT efficacy against hospitalization
    nat <- calculate_nat(variables = variables, index = hosp)
    infection_efficacy <- nat_efficacy_infection_cpp(ab_titre = nat,
                                                     parameters = parameters, 
                                                     day = day - 1L) # 0-based index
    severe_efficacy <- nat_efficacy_severe_cpp(ab_titre = nat,
                                               ef_infection = infection_efficacy,
                                               parameters = parameters,
                                               day = day - 1L) # 0-based index

    hosp$sample(prob_hosp * severe_efficacy)

    if (hosp$size() > availability) {
        hosp$choose(k = availability)
        availability <- 0
    } else {
        availability <- availability - hosp$size()
    }

    not_hosp <- target$set_difference(hosp)

    if (hosp$size() > 0) {
        events$enter_hospital_waiting$schedule(target = hosp, delay = 0)
        events$enter_hospital$schedule(target = hosp, delay = parameters$enter_hosp_delay)
        events$leave_hospital$schedule(target = hosp, delay = parameters$enter_hosp_delay +
                                                              parameters$dur_hospital)
    }

    disc_ages_hosp <- variables$discrete_age$get_values(hosp)
    disc_ages_not_hosp <- variables$discrete_age$get_values(not_hosp)

    prob_death_hosp <- rr_hosp * (1 - prob_recover[disc_ages_hosp])
    prob_death_not_hosp <- (1 - prob_recover[disc_ages_not_hosp])

    to_recover_hosp <- hosp$copy()
    to_recover_hosp$sample(1 - prob_death_hosp)
    to_dead_hosp <- hosp$set_difference(to_recover_hosp)

    to_recover_not_hosp <- not_hosp$copy()
    to_recover_not_hosp$sample(1 - prob_death_not_hosp)
    to_dead_not_hosp <- not_hosp$set_difference(to_recover_not_hosp)

    if (to_recover_hosp$size() > 0) {
        events$recovery$schedule(target = to_recover_hosp,
                                 delay = S_delay(n = to_recover_hosp$size()))
    }

    if (to_dead_hosp$size() > 0) {
        events$death$schedule(target = to_dead_hosp,
                              delay = D_delay(n = to_dead_hosp$size()))
    }

    if (to_recover_not_hosp$size() > 0) {
        events$recovery$schedule(target = to_recover_not_hosp,
                                 delay = S_delay(n = to_recover_not_hosp$size()))
    }

    if (to_dead_not_hosp$size() > 0) {
        events$death$schedule(target = to_dead_not_hosp,
                              delay = D_delay(n = to_dead_not_hosp$size()))
    }

    return(availability)
}


#' @title Determines severity of infection
#' @description Determines disease severity and outcome
#' @param events List of events in the model
#' @param variables List of variables in the model
#' @param parameters The model parameters
#' @param dt Size of time step
#' @param shift Schedule future events after minimum number of time steps delay
#' @export
create_infection_scheduler_listener <- function(events, variables, parameters, dt, shift = 0L) {

    # Make function to calculate population NAT

    calculate_nat <- make_calculate_nat(variables = variables)

    # Build delay functions

    S_delay_IAsymp <- make_const(parameters$dur_IAsymp, dt = dt, shift = shift)
    D_delay_IAsymp <- make_const(parameters$dur_IAsymp, dt = dt, shift = shift)

    S_delay_IMild <- make_const(parameters$dur_IMild , dt = dt, shift = shift)
    D_delay_IMild <- make_const(parameters$dur_IMild , dt = dt, shift = shift)

    S_delay_IModerate <- make_const(parameters$dur_IModerate, dt = dt, shift = shift)
    D_delay_IModerate <- make_const(parameters$dur_IModerate, dt = dt, shift = shift)

    S_delay_ISevere <- make_const(parameters$dur_ISevere, dt = dt, shift = shift)
    D_delay_ISevere <- make_const(parameters$dur_ISevere, dt = dt, shift = shift)

    S_delay_ICritical <- make_const(parameters$dur_ICritical, dt = dt, shift = shift)
    D_delay_ICritical <- make_const(parameters$dur_ICritical, dt = dt, shift = shift)

    return(
        function(timestep, target) {

            # Current day
            day <- ceiling(timestep * dt)

            prob_IAsymp <- parameters$prob_IAsymp
            prob_IMild <- parameters$prob_IMild
            prob_IModerate <- parameters$prob_IModerate
            prob_ISevere <- parameters$prob_ISevere
            prob_ICritical <- parameters$prob_ICritical

            # Determine disease severity
            disc_ages <- variables$discrete_age$get_values(target)
            prob <- prob_IAsymp[disc_ages]
            to_asymp <- target$copy()
            to_asymp$sample(prob)
            target <- target$set_difference(to_asymp)

            disc_ages <- variables$discrete_age$get_values(target)
            prob <- prob_IMild[disc_ages] / (prob_IMild[disc_ages] + prob_IModerate[disc_ages] +
                                             prob_ISevere[disc_ages] + prob_ICritical[disc_ages])
            to_mild <- target$copy()
            to_mild$sample(prob)
            target <- target$set_difference(to_mild)

            disc_ages <- variables$discrete_age$get_values(target)
            prob <- prob_IModerate[disc_ages] / (prob_IModerate[disc_ages] +
                                                 prob_ISevere[disc_ages] +
                                                 prob_ICritical[disc_ages])
            to_moderate <- target$copy()
            to_moderate$sample(prob)
            target <- target$set_difference(to_moderate)

            disc_ages <- variables$discrete_age$get_values(target)
            prob <- prob_ISevere[disc_ages] / (prob_ISevere[disc_ages] +
                                               prob_ICritical[disc_ages])
            to_severe <- target$copy()
            to_severe$sample(prob)
            to_critical <- target$set_difference(to_severe)

            # Determine available hospital beds

            waiting_for_hospital <- variables$hospitalization_states$get_index_of("Pre_Hospital")
            num_waiting_for_hosptial <- waiting_for_hospital$size()
            max_waiting_for_hosptial <- as.integer((parameters$enter_hosp_delay /
                                                    parameters$dur_hospital) *
                                                    parameters$beds_total)
            availability <- max_waiting_for_hosptial - num_waiting_for_hosptial

            # Determine outcomes in descending order of severity, thereby prioritizing availability

            availability <- determine_outcome(prob_recover = parameters$prob_ICritical_to_S,
                                              prob_hosp = parameters$prob_hosp_ICritical,
                                              rr_hosp = parameters$rr_hosp_ICritical,
                                              target = to_critical,
                                              availability = availability,
                                              S_delay = S_delay_ICritical,
                                              D_delay = D_delay_ICritical,
                                              events = events,
                                              parameters = parameters,
                                              calculate_nat = calculate_nat,
                                              day = day)

            availability <- determine_outcome(prob_recover = parameters$prob_ISevere_to_S,
                                              prob_hosp = parameters$prob_hosp_ISevere,
                                              rr_hosp = parameters$rr_hosp_ISevere,
                                              target = to_severe,
                                              availability = availability,
                                              S_delay = S_delay_ISevere,
                                              D_delay = D_delay_ISevere,
                                              events = events,
                                              parameters = parameters,
                                              calculate_nat = calculate_nat,
                                              day = day)

            availability <- determine_outcome(prob_recover = parameters$prob_IModerate_to_S,
                                              prob_hosp = parameters$prob_hosp_IModerate,
                                              rr_hosp = parameters$rr_hosp_IModerate,
                                              target = to_moderate,
                                              availability = availability,
                                              S_delay = S_delay_IModerate,
                                              D_delay = D_delay_IModerate,
                                              events = events,
                                              parameters = parameters,
                                              calculate_nat = calculate_nat,
                                              day = day)

            availability <- determine_outcome(prob_recover = parameters$prob_IMild_to_S,
                                              prob_hosp = parameters$prob_hosp_IMild,
                                              rr_hosp = parameters$rr_hosp_IMild,
                                              target = to_mild,
                                              availability = availability,
                                              S_delay = S_delay_IMild,
                                              D_delay = D_delay_IMild,
                                              events = events,
                                              parameters = parameters,
                                              calculate_nat = calculate_nat,
                                              day = day)

            availability <- determine_outcome(prob_recover = parameters$prob_IAsymp_to_S,
                                              prob_hosp = parameters$prob_hosp_IAsymp,
                                              rr_hosp = parameters$rr_hosp_IAsymp,
                                              target = to_asymp,
                                              availability = availability,
                                              S_delay = S_delay_IAsymp,
                                              D_delay = D_delay_IAsymp,
                                              events = events,
                                              parameters = parameters,
                                              calculate_nat = calculate_nat,
                                              day = day)

        }

    )

}


#' @title Infectiousness start time listener
#' @description Updates variable to record start time of an infectious period
#' @param infectiousness_start_time The start time variable
#' @export
create_infectiousness_start_time_update_listener <- function(infectiousness_start_time) {
    function(timestep, target) {
        infectiousness_start_time$queue_update(value = timestep, index = target)
    }
}


#' @title State update listener
#' @description Updates health states
#' @param states A \code{\link[individual]{CategoricalVariable}} object
#' @param destination The new health state
#' @export
create_state_update_listener <- function(states, destination) {
    function(timestep, target) {
        states$queue_update(value = destination, index = target)
    }
}


#' @title Hospitalization state update listener
#' @description Updates Hospitalization states
#' @param states A \code{\link[individual]{CategoricalVariable}} object
#' @param destination The new Hospitalization state
#' @export
create_hospitalization_state_update_listener <- function(hospitalization_states, destination) {
    function(timestep, target) {
        hospitalization_states$queue_update(value = destination, index = target)
    }
}


#' @title Attach listener functions
#' @description Attach listerner functions that respond to events
#' @param variables The list of model variables
#' @param events The list of [individual::TargetedEvent], the output of [create_events]
#' @param parameters The list of model parameters
#' @param dt Size of the time step
#' @param shift Size of delay shift
#' @export
attach_event_listeners <- function(variables, events, parameters, dt, shift = 0L) {

    # Infection ----------

    events$infection$add_listener(
        create_state_update_listener(
            variables$states,
            "I"
        )
    )

    events$infection$add_listener(
        create_infection_scheduler_listener(
            events,
            variables,
            parameters,
            dt = dt,
            shift = shift
        )
    )

    events$infection$add_listener(
        create_infectiousness_start_time_update_listener(
            variables$infectiousness_start_time
        )
    )

    # Recovery ----------
    events$recovery$add_listener(
        create_state_update_listener(
            variables$states,
            "S"
        )
    )

    # Death ----------
    events$death$add_listener(
        create_state_update_listener(
            variables$states,
            "D"
        )
    )

    # Enter hospital waiting ----------

    events$enter_hospital_waiting$add_listener(
        create_hospitalization_state_update_listener(
            variables$hospitalization_states,
            "Pre_Hospital"
        )
    )

    # Enter hospital ----------

    events$enter_hospital$add_listener(
        create_hospitalization_state_update_listener(
            variables$hospitalization_states,
            "In_Hospital"
        )
    )

    # Leave hospital ----------

    events$leave_hospital$add_listener(
        create_hospitalization_state_update_listener(
            variables$hospitalization_states,
            "Not_In_Hospital"
        )
    )

}


#' @title Schedule events for individuals at initialisation
#' @description Event initialization
#' @param parameters model parameters
#' @param events a list of [individual::TargetedEvent], the output of [create_events]
#' @param variables a list of model variables, the output of [create_variables]
#' @param dt size of the time step
#' @export
setup_events <- function(parameters, events, variables, dt) {

    # I
    bset_I <- variables$states$get_index_of("I")
    if (bset_I$size() > 0) {
        init_fn <- create_infection_scheduler_listener(
            events,
            variables,
            parameters,
            dt = dt
        )
        init_fn(timestep = 1, target = bset_I)
    }

}
