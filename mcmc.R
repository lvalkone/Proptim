# Don't run directly, called by source() from model_*.R files
mcmc_wrapper <- function(seed, code, use_all_data, data_subs, data_subsC, 
                         data_nonsC, data_nonsC_never, id_collector_subs,
                         id_collector_nonsubs, monitors, niter, nburnin, thin) {
  
  library("nimble")
  library("nimbleHMC")
  library("magrittr")
  library("dplyr")
  
  set.seed(seed)
  
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
      use_all_data = use_all_data,
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
    alpha1 = rnorm(1, 0.35, 0.01),
    alpha2 = rnorm(1, 0.45, 0.01),
    alpha3 = rnorm(1, -0.30, 0.01),
    beta_age = c(0.0, rnorm(N_A - 1, c(-0.15, -0.30, -0.45), 0.01)),
    beta_const = rnorm(1, 2.8, 0.01),
    beta_gender = c(0.0, rnorm(N_G - 1, 0.01, 0.01)),
    beta_location = c(0.0, rnorm(N_L - 1, -0.02, 0.01)),
    kappa = rnorm(1, 0.75, 0.01),
    u_nonC_never_raw = rnorm(N_nonC_never, 0.0, 0.01),
    u_sub_raw = rnorm(N_sub, 0.0, 0.01),
    u_tau = rnorm(1, 0.1, 0.01)
  ))

  # Model definition
  q_model <- nimbleModel(
    code = code,
    constants = q_const,
    data = q_data,
    inits = q_init,
    buildDerivs = TRUE
  )

  # Use default HMC configuration
  # Don't print due to large number of parameters
  q_conf <- configureHMC(q_model, print = FALSE)

  # Set which parameters to monitor
  q_conf$setMonitors(monitors)

  # Generates the MCMC functions to be used for sampling
  q_mcmc <- buildMCMC(q_conf)

  # Compile the model
  q_Cmodel <- compileNimble(q_model)

  # Compile the functions used for MCMC
  q_Cmcmc <- compileNimble(q_mcmc)

  # Get the posterior samples
  out <- runMCMC(
    q_Cmcmc,
    niter = niter,
    nburnin = nburnin,
    thin = thin,
    setSeed = seed
  )

  out
}

# Run MCMC with parallel chains
nc <- min(detectCores(), 12)
clust <- makeCluster(nc)
samples <- parLapply(
  cl = clust, 
  X = seq_len(nc),
  fun = mcmc_wrapper,
  code = q_modelcode,
  use_all_data = TRUE,
  data_subs = data_subs,
  data_subsC = data_subsC,
  data_nonsC = data_nonsC, 
  data_nonsC_never = data_nonsC_never,
  id_collector_subs = id_collector_subs,
  id_collector_nonsubs = id_collector_nonsubs,
  monitors = monitors,
  niter = 22000,
  nburnin = 2000,
  thin = 50
)
stopCluster(clust)

saveRDS(samples, file = paste0("posterior_samples", suffix, ".rds"))
