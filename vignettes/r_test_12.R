library(GA)

time_period <- 100

run <- function(beta_set)
{
  target <- rep(0.5, time_period)
  1 / (sum((beta_set - target) ^ 2))
}

GA <- ga(type = "real-valued",
         fitness = function(x) run(x),
         lower = rep(0.0, time_period), upper = rep(1.0, time_period),
         popSize = 50, maxiter = 1000, run = 100, parallel = 5)
summary(GA)
