# Notes

1. If we need to include onward infectiousness as something dependent upon ab titre
for breakthrough infections, it will need to interact with this: 
```
# Group infection by age
ages <- variables$discrete_age$get_values(infectious)
inf_ages <- tab_bins(a = ages, nbins = parameters$number_of_age_groups)
```
If it can be turned into a vector of values of relative infectiousness the same
length as `ages` then we just need to make a modified `tab_bins` that handles that
relative multiplicative contribution of each individual.

2. One thing to note is that when a large amount of doses are available from the
start of the simulation there may be wasted doses during the lag between people
getting their first dose, and when persons can get (become eligible for) their second dose.
The prioritization step will only increment when both are fulfilled.

3. one way to include additional variation in inter-dose times is to 
include variance in `dose_period` used to calculate thresholds in 
`threshold <- timestep - as.integer(parameters$dose_period[dose]/dt)` from
function `distribution_vaccines.R/get_current_eligible_from_coverage`.

4. including multiple vaccine types could be as simple as just finding out what
proportion of doses belong to each type per time step and draw the ab titre
from the appropriate parameter.
