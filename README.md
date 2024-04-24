## Overview

Code files and data set related to preprint "Price Optimization Combining Conjoint Data and Purchase History: A Causal Modeling Approach".


## Running instructions

For running the workflow with R software:

1. `simulation.R`
    - Creates the population using `finns_agl.csv` (Statistics Finland (2021))
    - Simulates the purchase history and conjoint studies
    - Saves the simulated data and the true expected gross profits into `.RData` files
2. `dosearch_transportability.R`
    - Searches the identifiability of the causal effect of price on the purchase choice
3. `model.R`
    - Fits a Bayesian model for the purchase choice
    - Saves the obtained posterior samples and statistics into `.RData` files
4. `output.R`
    - Generates posterior predictive samples for calculating the estimated expected gross profit
    - Compares the estimated and the true expected gross profits

### Additional files

- `model_misspecified_Q.R`, `model_misspecified_kappa.R`, `model_misspecified_prds.R`
    - Model codes for different misspecification scenarios 
- `mcmc.R`
    - Code for running the model (called from the `.R` files above)

## References

Statistics Finland (2021). Population 31.12. by Information, Urban-rural classification,
Age, Area, Year and Sex. https://pxdata.stat.fi/PxWeb/pxweb/en/StatFin/. The
material was downloaded from Statistics Finland’s interface service on 2021-09-27 with
the licence CC BY 4.0.
