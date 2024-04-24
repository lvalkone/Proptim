library("nimble")
library("nimbleHMC")
library("tibble")
library("dplyr")
library("parallel")

load(file = "cabc_data.RData") # Population data, purchase history and conjoint data

set.seed(0)

# In this model, the conjoint effect is misspecified as multiplicative instead of additive

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
    logit(p_subC[i]) <- alpha1 * (T_subC[Id_subC[i]] - kappa * X_subC[i]) +
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
    logit(p_nonC[i]) <- alpha1 * (T_nonC[Id_nonC[i]] - kappa * X_nonC[i]) +
      alpha2 * log1p(Prds_nonC[i]) + alpha3 * (Prds_nonC[i] == 0)
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
  for (i in 1:(M_nonC_never)) {
    logit(p_nonC_never[i]) <- alpha1 * (T_nonC_never[Id_nonC_never[i]] - kappa * X_nonC_never[i]) +
      alpha2 * log1p(Prds_nonC_never[i]) + alpha3 * (Prds_nonC_never[i] == 0)
      Y_nonC_never[i] ~ dbin(prob = p_nonC_never[i], size = 1)
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
  kappa ~ dgamma(5, 5)
  alpha1 ~ dnorm(0, sd = sd_init)
  alpha2 ~ dnorm(0, sd = sd_init)
  alpha3 ~ dnorm(0, sd = sd_init)
  u_tau ~ dgamma(shape = shape, scale = scale)
})

# Model constants
q_const <- c(
  with(data_subs, list(
    N_sub = subjects %>% nrow(),
    M_sub = tasks %>% nrow(),
    A_sub = as.numeric(subjects$a),
    G_sub = as.numeric(subjects$g),
    L_sub = as.numeric(subjects$l),
    Prds_sub = as.numeric(tasks$prds),
    Id_sub = as.numeric(tasks$uid),
    N_A = subjects %>% select(a) %>% n_distinct(na.rm = TRUE),
    N_G = subjects %>% select(g) %>% n_distinct(na.rm = TRUE),
    N_L = subjects %>% select(l) %>% n_distinct(na.rm = TRUE)
  )),
  with(data_subsC, list(
    id_collector_subs = as.numeric(id_collector_subs),
    N_subC = subjects %>% nrow(),
    M_subC = tasks %>% nrow(),
    Prds_subC = as.numeric(tasks$prds),
    Id_subC = as.numeric(tasks$uid)
  )),
  with(data_nonsC, list(
    id_collector_nonsubs = as.numeric(id_collector_nonsubs),
    N_nonC = subjects %>% nrow(),
    M_nonC = tasks %>% nrow(),
    Prds_nonC = as.numeric(tasks$prds),
    Id_nonC = as.numeric(tasks$uid)
  )),
  with(data_nonsC_never, list(
    N_nonC_never = subjects %>% nrow(),
    M_nonC_never = tasks %>% nrow(),
    A_nonC_never = as.numeric(subjects$a),
    G_nonC_never = as.numeric(subjects$g),
    L_nonC_never = as.numeric(subjects$l),
    Prds_nonC_never = as.numeric(tasks$prds),
    Id_nonC_never = as.numeric(tasks$uid)
  )),
  list(
    sd_init = 0.5,
    shape = 2,
    scale = 0.2
  )
)

suffix <- "_misspecified_kappa"
monitors <- c(
  "beta_const", "beta_age", "beta_gender", "beta_location",
  "kappa", "u_sub", "u_nonC_never", "u_tau", "alpha1", "alpha2", "alpha3"
)
source("MCMC.R")
