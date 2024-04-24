library("nimble")
library("nimbleHMC")
library("tibble")
library("dplyr")
library("parallel")

load(file = "cabc_data.RData") # Population data, purchase history and conjoint data

set.seed(0)

# In this model, the purchase probability is misspecified by omitting
# the effects of the earlier subscription periods

# Model code
q_modelcode <- nimbleCode({
  # Likelihoods
  # Subscription data
  # Reference price model
  for (i in 1:N_sub) {
    u_sub_raw[i] ~ dnorm(0, 1)
    u_sub[i] <- u_tau * u_sub_raw[i]
    T_sub[i] <- exp(
      beta_const +
        beta_age[A_sub[i]] +
        beta_gender[G_sub[i]] +
        beta_location[L_sub[i]] +
        u_sub[i]
    )
  }
  # Purchase choices
  for (i in 1:M_sub) {
    logit(p_sub[i]) <- alpha1 * (T_sub[Id_sub[i]] - X_sub[i])
    Y_sub[i] ~ dbin(prob = p_sub[i], size = 1)
  }
  if (use_all_data) {
    # Conjoint (Current subscribers)
    # Reference price model
    for (i in 1:N_subC) {
      T_subC[i] <- T_sub[id_collector_subs[i]]
    }
    # Purchase choices
    for (i in 1:M_subC) {
      logit(p_subC[i]) <- alpha1 * (T_subC[Id_subC[i]] + kappa - X_subC[i])
      Y_subC[i] ~ dbin(prob = p_subC[i], size = 1)
    }
    # Conjoint (Earlier subscribers)
    # Reference price model
    for (i in 1:N_nonC) {
      T_nonC[i] <- T_sub[id_collector_nonsubs[i]]
    }
    # Purchase choices
    for (i in 1:M_nonC) {
      logit(p_nonC[i]) <- alpha1 * (T_nonC[Id_nonC[i]] + kappa - X_nonC[i])
      Y_nonC[i] ~ dbin(prob = p_nonC[i], size = 1)
    }
    # Conjoint (Non-subscribers; never been customers)
    # Reference price model
    for (i in 1:N_nonC_never) {
      u_nonC_never_raw[i] ~ dnorm(0, 1)
      u_nonC_never[i] <- u_tau * u_nonC_never_raw[i]
      T_nonC_never[i] <- exp(
        beta_const +
          beta_age[A_nonC_never[i]] +
          beta_gender[G_nonC_never[i]] +
          beta_location[L_nonC_never[i]] +
          u_nonC_never[i]
      )
    }
    # Purchase choices
    for (i in 1:M_nonC_never) {
      logit(p_nonC_never[i]) <- alpha1 * (T_nonC_never[Id_nonC_never[i]] + kappa - X_nonC_never[i])
      Y_nonC_never[i] ~ dbin(prob = p_nonC_never[i], size = 1)
    }
  }
  # Priors
  beta_const ~ dnorm(0, sd = sd_init)
  beta_age[1] <- 0
  for (i in 2:N_A) {
    beta_age[i] ~ dnorm(0, sd = sd_init)
  }
  beta_gender[1] <- 0
  beta_gender[2] ~ dnorm(0, sd = sd_init)
  beta_location[1] <- 0
  beta_location[2] ~ dnorm(0, sd = sd_init)
  if (use_all_data) {
    kappa ~ dnorm(0, sd = sd_init)
  }
  alpha1 ~ dnorm(0, sd = sd_init)
  u_tau ~ dgamma(shape = shape, scale = scale)
})

suffix <- "_misspecified_prds"
monitors <- c(
  "beta_const", "beta_age", "beta_gender", "beta_location", 
  "kappa", "u_sub", "u_nonC_never", "u_tau", "alpha1"
)
source("MCMC.R")
