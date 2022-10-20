# **safir**

<!-- badges: start -->
[![R build
status](https://github.com/mrc-ide/safir/workflows/R-CMD-check/badge.svg)](https://github.com/mrc-ide/safir/actions)
[![CodeFactor](https://www.codefactor.io/repository/github/mrc-ide/safir/badge)](https://www.codefactor.io/repository/github/mrc-ide/safir)
[![codecov.io](https://codecov.io/github/mrc-ide/safir/coverage.svg?branch=main)](https://codecov.io/github/mrc-ide/safir?branch=main)
<!-- badges: end -->

**safir**: **s**quire **a**nd **f**riends **i**ndividual **r**ewrite
(and maintains the squire naming
[theme](https://en.wikipedia.org/wiki/Knights_of_the_Round_Table#Safir)).
It uses the [`{individual}`](https://github.com/mrc-ide/individual)
package to specify and run the simulation.

**safir** can run individual based stochastic versions of the **squire** and
**nimue** models, which recover trajectories from those aggregated population
models (see package vignettes). It also implements a model for vaccination
with an arbitrary number of doses and antibody titre dynamics which affect
protective efficacy as well as efficacy against severe disease outcomes.

Because **safir** implements multiple models as individual based simulations,
there is a large amount of code, but the central disease progression is shared
across the models within **safir**. Each vignette describes how the simulation
functions and where the relevant code can be found.

## Installation

``` r
install_github('mrc-ide/safir')
library(safir)
```

## Documentation

[`{safir}`](https://github.com/mrc-ide/safir) is documented on a
[dedicated website](https://mrc-ide.github.io/safir).

This includes the following vignettes:

-   **`Squire Validation Run`**: comparison of
    [`{safir}`](https://github.com/mrc-ide/safir) to
    [`{squire}`](https://github.com/mrc-ide/squire).
-   **`Nimue Validation Run`**: comparison of
    [`{safir}`](https://github.com/mrc-ide/safir) to
    [`{nimue}`](https://github.com/mrc-ide/nimue).
-   **`Vaccine model with multiple doses`**: runs the model for
    explicit tracking of antibody titre for a 2 dose vaccine.
-   **`Natural Immunity`**: demonstrates how to modify the vaccine model
    to also include mechanistic incorporation of natural immunity
    from infection, as opposed to a generic R (recovered) class.
-   **`Modeling of antibody titre and immunity`**: describes how antibody
    titre and protection against infection and severe disease is modeled.
-   **`Differential modeling of infection and vaccine derived NATs`**: is a short
    description of how to model differential boost/decay of vaccine and infection
    derived neutralizing antibody titres.
    
## License

MIT © Imperial College of Science, Technology and Medicine
