# **safir3**

**safir** runs individual based stochastic models. It is based on the
[`{individual}`](https://github.com/mrc-ide/individual) package. It can be
used to model vaccination, with an arbitrary number of doses and antibody
titre dynamics which affect protective efficacy as well as efficacy
against severe disease outcomes.

**safir3** adds individual infectiousness, for modelling the impact of
antiviral treatments, a simplified compartmental system, and a few other features.

## Documentation

The original [`{safir}`](https://github.com/mrc-ide/safir) is documented on a
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
