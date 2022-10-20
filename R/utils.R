#' @noRd
get_probabilties <- function(prob, ages, day) {
  if (inherits(prob, "matrix")) {
    prob[cbind(ages, day)]
  } else {
    prob[ages]
  }
}

#' @noRd
get_size_bset <- function(bsets) {
  vapply(X = bsets, FUN = function(x) {x$size()}, FUN.VALUE = numeric(1), USE.NAMES = FALSE)
}

#' @noRd
check_probabilities <- function(prob, parameters) {
  if (inherits(prob, "matrix")) {
    stopifnot(nrow(prob) == parameters$number_of_age_groups)
    stopifnot(ncol(prob) == parameters$time_period)
  } else {
    stopifnot(length(prob) == parameters$number_of_age_groups)
  }
  stopifnot(is.finite(prob))
}

`%||%` <- function(a, b) { # nolint
  if (is.null(a)) b else a
}

#' @title Contact matrix at current timestep
#' @description Currently only supports one contact matrix
#' @param parameters object returned from [get_parameters]
#' @return 2D contact matrix for the current timestep
#' @noRd
get_contact_matrix <- function(parameters) {
  parameters$mix_mat_set[1, , ]
}

#' @title Produce piecewise linear Rt
#' @description Make a piecewise linear interpolating vector of Rt, for each
#' day within the simulated interval. The assumed start of simulation
#' is the first date provided in `dates`, and the end of simulation is the last
#' date of `dates` unless `max_date` is provided. Interpolation is done with
#' inclusive endpoints.
#' @param dates a vector of [Date]s
#' @param rt a vector of reproductive values at those dates
#' @param max_date the maximum date of simulation (optional)
#' @return a [list] with named elements `Rt` and `Rt_tt`, which can be provided
#' to the [safir3::get_parameters] function's arguments `R0` and `tt_R0`.
#' @importFrom stats approx
#' @export
interpolate_rt <- function(dates, rt, max_date = NULL) {
  stopifnot(inherits(dates, "Date"))
  stopifnot(length(dates) > 1)
  stopifnot(all(as.integer(diff(dates)) > 0))
  stopifnot(is.finite(rt))
  stopifnot(all(rt > 0))
  stopifnot(length(dates) == length(rt))

  if (!is.null(max_date)) {
    stopifnot(inherits(max_date, "Date"))
    stopifnot(max_date > dates[length(dates)])
    dates <- c(dates, max_date)
    rt <- c(rt, rt[length(rt)])
  }

  time_interval <- vapply(X = 2:length(dates), FUN = function(d){
    as.integer(difftime(dates[d], dates[1] - 1))
  }, FUN.VALUE = integer(1))
  time_interval <- c(0, time_interval)

  total_interval <- as.integer(difftime(dates[length(dates)], dates[1] - 1))
  Rt_tt <- 1:total_interval
  Rt <- rep(NaN, total_interval)

  time_interval_list <- lapply(X = seq_len(length(time_interval) - 1), FUN = function(v){
    (time_interval[v] + 1):time_interval[v + 1]
  })

  for (i in seq_len(length(dates) - 1)) {

    if (rt[i] == rt[i + 1]) {
      Rt[time_interval_list[[i]]] <- rt[i]
    } else {
      Rt[time_interval_list[[i]]] <- approx(
        x = c(time_interval_list[[i]][1], time_interval_list[[i]][length(time_interval_list[[i]])]),
        y = c(rt[i], rt[i + 1]),
        xout = time_interval_list[[i]]
      )$y
    }

  }

  return(list(
    Rt = Rt, Rt_tt = Rt_tt
  ))
}

#' @title Make Erlang waiting time distribution
#'
#' @description Random draws from Erlang distribution
#' scaled by time step size.
#' @param mu Mean duration
#' @param dt Size of time step
#' @param shape Shape parameter of Erlang distribution
#' @param shift number of time steps to add to sampled value
#' @importFrom stats rgamma
#' @export
make_rerlang <- function(mu, dt, shape = 2, shift = 0L) {
  assert_pos(mu, zero_allowed = FALSE)
  assert_pos(dt, zero_allowed = FALSE)
  assert_pos(shape, zero_allowed = FALSE)
  r <- shape / mu
  function(n) {
    floor(rgamma(n = n, shape = shape, rate = r) / dt) + shift # possibility of 0 delay
  }
}


#' @title Make discretized exponential waiting time distribution
#' @description Make geometric approximation to
#' continuous time exponential distribution.
#' @param mu Mean duration
#' @param dt Size of time step
#' @param shift number of time steps to add to sampled value
#' @importFrom stats rgeom pexp
make_rexp <- function(mu, dt, shift = 0L) {
  assert_pos(mu, zero_allowed = FALSE)
  assert_pos(dt, zero_allowed = FALSE)
  r <- 1 / mu
  p <- pexp(q = r * dt)
  function(n) {
    rgeom(n = n, prob = p) + shift
  }
}


#' @title Make constant waiting time function
#' @description Return constant waiting time
#' @param me Fixed duration
#' @param dt Size of time step
#' @param shift number of time steps to add to sampled value
make_const <- function(mu, dt, shift = 0L) {
  assert_pos(mu, zero_allowed = FALSE)
  assert_pos(dt, zero_allowed = FALSE)
  function(n) {
    mu + shift
  }
}

#' @noRd
remove_non_numerics <- function(l) {
  clean <- list()
  for (key in names(l)) {
    if (storage.mode(l[[key]]) %in% c("integer", "double")) {
      clean[[key]] <- l[[key]]
    }
  }
  clean
}


#' @noRd
vcapply <- function(X, FUN, ...) vapply(X, FUN, ..., character(1))


#'@noRd
rep_list <- function(x, times) {
  setNames(lapply(x, function(x) {rep(list(x), times)}), names(x))
}
