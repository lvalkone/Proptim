
library(nimble)
library(tibble)
library(dplyr)
library(Rcpp)
library(coda)

rm(list=ls())
load(file = "cabcData.RData") # Population data, purchase history and conjoint data
set.seed(2492021)


########################################################################################
#   FIT THE MODEL
########################################################################################

# Indicator: FALSE means we use only sales data
use_all_data <- TRUE


# Model code

q_modelcode <- nimbleCode({
  
  # Likelihoods
  
  # Subscription data
  
  # Reference price model
  for (i in 1:N_sub) {
    
    u_sub[i] ~ dnorm(0, sd = u_tau)
    
    T_sub[i] <- exp(beta_const +
                      beta_age[A_sub[i]] +
                      beta_gender[G_sub[i]] +
                      beta_location[L_sub[i]] + 
                      u_sub[i])
    
  }
  
  # Purchase choices
  for (i in 1:(M_sub)) {
    
    logit(p_sub[i]) <- alpha1*(T_sub[Id_sub[i]] - X_sub[i]) +
      alpha2*log(Prds_sub[i]+1) + alpha3*(Prds_sub[i]==0)
    
    
    Y_sub[i] ~ dbin(prob = p_sub[i], size = 1)
    
  }
  

  if(use_all_data == TRUE) {
    
    # Conjoint (Current subscribers)
    
    # Reference price model
    for (i in 1:N_subC) {
      
      T_subC[i] <- T_sub[id_collector_subs[i]]
      
    }
    
    # Purchase choices
    for (i in 1:M_subC) {
      
      logit(p_subC[i]) <- alpha1*(T_subC[Id_subC[i]] + kappa - X_subC[i]) +
        alpha2*log(Prds_subC[i]+1) + alpha3*(Prds_subC[i]==0)
      
      Y_subC[i] ~ dbin(prob = p_subC[i], size = 1)
      
    }
    
    
    # Conjoint (Earlier subscribers)
    
    # Reference price model
    for (i in 1:N_nonC) {
      
      T_nonC[i] <- T_sub[id_collector_nonsubs[i]]
      
    }
    
    # Purchase choices
    for (i in 1:M_nonC) {
      
      logit(p_nonC[i]) <- alpha1*(T_nonC[Id_nonC[i]] + kappa - X_nonC[i]) +
        alpha2*log(Prds_nonC[i]+1) + alpha3*(Prds_nonC[i]==0)
      
      Y_nonC[i] ~ dbin(prob = p_nonC[i], size = 1)
      
    }
    
    
    # Conjoint (Non-subscribers; never been customers)
    
    # Reference price model
    for (i in 1:N_nonC_never) {
      
      u_nonC_never[i] ~ dnorm(0, sd = u_tau)
      
      T_nonC_never[i] <- exp(beta_const +
                               beta_age[A_nonC_never[i]] +
                               beta_gender[G_nonC_never[i]] +
                               beta_location[L_nonC_never[i]] + 
                               u_nonC_never[i])
      
    }
    
    # Purchase choices
    for (i in 1:(M_nonC_never)) {
      
      logit(p_nonC_never[i]) <- alpha1*(T_nonC_never[Id_nonC_never[i]] + kappa - X_nonC_never[i]) +
        alpha2*log(Prds_nonC_never[i]+1) + alpha3*(Prds_nonC_never[i]==0)

      Y_nonC_never[i] ~ dbin(prob = p_nonC_never[i], size = 1)
      
    }

  }
  
  # Priors
  
  # betas
  beta_const ~ dnorm(0, sd = sd_init)
  
  beta_age[1] <- 0
  for (i in 2:N_A) {
    beta_age[i] ~ dnorm(0, sd = sd_init)
  }
  
  beta_gender[1] <- 0
  beta_gender[2] ~ dnorm(0, sd = sd_init)
  
  beta_location[1] <- 0
  beta_location[2] ~ dnorm(0, sd = sd_init)
  
  # kappa
  if(use_all_data==TRUE) {
    kappa ~ dnorm(0, sd = sd_init)
  }
  
  # alphas
  alpha1 ~ dnorm(0, sd = sd_init)
  alpha2 ~ dnorm(0, sd = sd_init)
  alpha3 ~ dnorm(0, sd = sd_init)
  
  # u_tau
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


# Model data

q_data <- c(
  # Subscription data
  with(data_subs, list(
    Y_sub = as.numeric(tasks$y),
    X_sub = as.numeric(tasks$x1)
  )),
  
  # Conjoint (subscribers data)
  with(data_subsC, list(
    Y_subC = as.numeric(tasks$y),
    X_subC = as.numeric(tasks$x1)
  )),
  
  # Conjoint (non-subscribers data; earlier)
  with(data_nonsC, list(
    Y_nonC = as.numeric(tasks$y),
    X_nonC = as.numeric(tasks$x1)
  )),
  
  # Conjoint (non-subscribers data; never)
  with(data_nonsC_never, list(
    Y_nonC_never = as.numeric(tasks$y),
    X_nonC_never = as.numeric(tasks$x1)
  ))
)


# Initial values

q_init <- with(q_const, list(
  beta_const = 1,
  beta_age = rep(0, N_A),
  beta_gender = rep(0, N_G),
  beta_location = rep(0, N_L),
  kappa = 0,
  u_tau = 0.1,
  u_sub = rep(0.0, N_sub),
  T_sub = rep(0.1, q_const$N_sub),
  T_subC = rep(0.1, q_const$N_subC),
  T_nonC = rep(0.1, q_const$N_nonC),
  alpha1 = 0,
  alpha2 = 0,
  alpha3 = 0,
  p_sub = rep(0.0, M_sub),
  p_subC = rep(0.0, M_subC),  
  p_nonC = rep(0.0, M_nonC),
  u_nonC_never = rep(0.0, N_nonC_never),
  T_nonC_never = rep(0.1, q_const$N_nonC_never),
  p_nonC_never = rep(0.0, M_nonC_never)
))


# Model definition
q_model <- nimbleModel(code = q_modelcode,
                       constants = q_const,
                       data = q_data,
                       inits = q_init)


# Use default MCMC configuration
q_conf <- configureMCMC(q_model)


# Which parameters to monitor
if(use_all_data == TRUE) {
  q_conf$addMonitors(c("beta_const", "beta_age", "beta_gender", "beta_location",
                       "kappa", "u_sub", "u_nonC_never", "u_tau", "alpha1", "alpha2", "alpha3"))
}else{ #using only sales data
  q_conf$addMonitors(c("beta_const", "beta_age", "beta_gender", "beta_location",
                       "u_sub", "u_tau", "alpha1", "alpha2", "alpha3"))
}


# Generates the MCMC functions to be used for sampling
q_mcmc <- buildMCMC(q_conf)


# Compile the model
q_Cmodel <- compileNimble(q_model)


# Compile the functions used for MCMC
q_Cmcmc <- compileNimble(q_mcmc, project = q_model)


# Run the model
niter <- 1000000
nburnin <- 500000
thin <- 500
q_Cmcmc$run(niter = niter, nburnin = nburnin, thin = thin)


# Get the posterior samples
q_samples <- as.list(q_Cmcmc$mvSamples)


# MCMC convergence diagnostics
model_output <- as.matrix(q_Cmcmc$mvSamples)
model_output <- mcmc(model_output)
plot(model_output[,c(1:13,length(colnames(model_output)))])


# Save estimates
if(use_all_data==TRUE){
  
  # Save the posterior estimates
  post_means <- apply(model_output[,c(1:13,ncol(model_output))], 2, mean)
  post_sds <- apply(model_output[,c(1:13,ncol(model_output))], 2, sd)
  
  # Save posterior samples as a list (excluding u_sub and u_nonC_never)
  posterior <- list(
    model_output[,1],
    model_output[,2],
    model_output[,3],
    model_output[,4:7],
    model_output[,8],
    model_output[,9:10],
    model_output[,11:12],
    model_output[,13],
    model_output[,ncol(model_output)]
  )
  
  names(posterior) <- c(
    "alpha1",
    "alpha2",
    "alpha3",
    "beta_age",
    "beta_const",
    "beta_gender",
    "beta_location",
    "kappa",
    "u_tau"
  )
  
  
  # Save posterior samples u_nonC_never as a list
  posterior_u_nonC_never <- list(
    model_output[,grepl("u_nonC_never", colnames(model_output)) %>% which()]
  )
  
  names(posterior_u_nonC_never) <- c(
    "u_nonC_never"
  )
  
  
  # Save posterior samples u_sub as a list
  posterior_u <- list(
    model_output[,grepl("u_sub", colnames(model_output)) %>% which()]
  )
  
  names(posterior_u) <- c(
    "u_sub"
  )
  
  # Save the data and posterior samples
  save(posterior, posterior_u, posterior_u_nonC_never, file = "post_samples.RData")
  save(post_means, post_sds, file = "posterior_statistics.RData")
  
}else { #using only sales data
  
  # Save the posterior estimates
  post_means <- apply(model_output[,c(1:12,ncol(model_output))], 2, mean)
  post_sds <- apply(model_output[,c(1:12,ncol(model_output))], 2, sd)
  
  # Save posterior samples as a list (excluding u_sub)
  posterior <- list(
    model_output[,1],
    model_output[,2],
    model_output[,3],
    model_output[,4:7],
    model_output[,8],
    model_output[,9:10],
    model_output[,11:12],
    model_output[,ncol(model_output)]
  )
  
  names(posterior) <- c(
    "alpha1",
    "alpha2",
    "alpha3",
    "beta_age",
    "beta_const",
    "beta_gender",
    "beta_location",
    "u_tau"
  )

  # Save posterior samples u_sub as a list
  posterior_u <- list(
    model_output[,grepl("u_sub", colnames(model_output)) %>% which()]
  )
  
  names(posterior_u) <- c(
    "u_sub"
  )

  # Save the data and posterior samples
  save(posterior, posterior_u, file = "post_samples_only_sales.RData")
  save(post_means, post_sds, file = "posterior_statistics_only_sales.RData")
  
}

