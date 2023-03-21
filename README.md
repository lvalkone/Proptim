## Overview

Code files and data set related to preprint "Price Optimization Combining Conjoint Data and Purchase History: A Causal Modeling Approach".

Preprint available here (will be added)


## Running instructions

For running the workflow with R software:

1. `data_simulation.R`
    - Creates the population using `finns_agl.csv`
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

