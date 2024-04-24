library("nimble")
library("nimbleHMC")
library("tibble")
library("dplyr")
library("parallel")

load(file = "cabc_data.RData") # Population data, purchase history and conjoint data

set.seed(0)

# In this model, Q is misspecified such that the effect of age, gender and location
# are not accounted for

# Model code
q_modelcode <- nimbleCode({
  # Likelihoods
  # Subscription data
  # Reference price model
  for (i in 1:N_sub) {
    u_sub_raw[i] ~ dnorm(0, 1)
    u_sub[i] <- u_tau * u_sub_raw[i]
    T_sub[i] <- exp(beta_const + u_sub[i])
  }
  # Purchase choices
  for (i in 1:M_sub) {
    logit(p_sub[i]) <- alpha1 * (T_sub[Id_sub[i]] - X_sub[i]) +
      alpha2 * log1p(Prds_sub[i]) + alpha3 * (Prds_sub[i] == 0)
    Y_sub[i] ~ dbin(prob = p_sub[i], size = 1)
  }
  # Conjoint (Current subscribers)
  # Reference price model
  for (i in 1:N_subC) {
    T_subC[i] <- T_sub[id_collector_subs[i]]
  }
  # Purchase choices
  for (i in 1:M_subC) {
    logit(p_subC[i]) <- alpha1 * (T_subC[Id_subC[i]] + kappa - X_subC[i]) +
      alpha2 * log1p(Prds_subC[i]) + alpha3 * (Prds_subC[i] == 0)
    Y_subC[i] ~ dbin(prob = p_subC[i], size = 1)
  }
  # Conjoint (Earlier subscribers)
  # Reference price model
  for (i in 1:N_nonC) {
    T_nonC[i] <- T_sub[id_collector_nonsubs[i]]
  }
  # Purchase choices
  for (i in 1:M_nonC) {
    logit(p_nonC[i]) <- alpha1 * (T_nonC[Id_nonC[i]] + kappa - X_nonC[i]) +
      alpha2 * log1p(Prds_nonC[i]) + alpha3 * (Prds_nonC[i] == 0)
    Y_nonC[i] ~ dbin(prob = p_nonC[i], size = 1)
  }
  # Conjoint (Non-subscribers; never been customers)
  # Reference price model
  for (i in 1:N_nonC_never) {
    u_nonC_never_raw[i] ~ dnorm(0, 1)
    u_nonC_never[i] <- u_tau * u_nonC_never_raw[i]
    T_nonC_never[i] <- exp(beta_const + u_nonC_never[i])
  }
  # Purchase choices
  for (i in 1:(M_nonC_never)) {
    logit(p_nonC_never[i]) <- alpha1 * (T_nonC_never[Id_nonC_never[i]] + kappa - X_nonC_never[i]) +
      alpha2 * log1p(Prds_nonC_never[i]) + alpha3 * (Prds_nonC_never[i] == 0)
      Y_nonC_never[i] ~ dbin(prob = p_nonC_never[i], size = 1)
  }
  # Priors
  beta_const ~ dnorm(0, sd = sd_init)
  kappa ~ dnorm(0, sd = sd_init)
  alpha1 ~ dnorm(0, sd = sd_init)
  alpha2 ~ dnorm(0, sd = sd_init)
  alpha3 ~ dnorm(0, sd = sd_init)
  u_tau ~ dgamma(shape = shape, scale = scale)
})

suffix <- "_misspecified_Q"
monitors <- c(
  "beta_const", "kappa", "u_sub", "u_nonC_never", 
  "u_tau", "alpha1", "alpha2", "alpha3"
)
source("MCMC.R")
