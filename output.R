library("dplyr")
library("magrittr")
library("ggplot2")
library("patchwork")

load("real_optimum.RData")
load("cabc_data.RData")
load("newdata.RData")

set.seed(0)

# indicator: FALSE means we use only sales data (only for the main model)
use_all_data <- TRUE

suffix <- ""

# If using correct model
samples <- if (use_all_data) {
  readRDS("posterior_samples.rds")
} else {
  readRDS("posterior_samples_only_sales.rds")
}

# Samples for misspecified models:

# Misspecified reference price
# samples <- readRDS("posterior_samples_misspecified_Q.rds")
# suffix <- "_misspecified_Q"

# Misspecified conjoint effect
# NOTE: Convergence issues with chains 5 and 7, so we drop these here
# Effective sample sizes are sufficient with the remaining chains
# samples <- readRDS("posterior_samples_misspecified_kappa.rds")
# samples <- samples[c(-5, -7)]
# suffix <- "_misspecified_kappa"

# Misspecified subscription periods
# NOTE: Convergence issues with chain 6, se we drop this
# Effective sample sizes are sufficient with the remaining chains
# samples <- readRDS("posterior_samples_misspecified_preds.rds")
# samples <- samples[-6]
# suffix <- "_misspecified_prds"


# Preparations for the t = 25 data ----------------------------------------


# New customers

# Initialize an index showing the data source of the observation
# 0 = completely new customers
# 1 = earlier customers
# 2 = conjoint of 'never been customers'
pop_H0$source <- 0

# 1) Those new customers in H0 who has been earlier subscribers
ind_sub <- intersect(pop_H0$uido, data_subs$subjects$uido)

# 2) Those new customers in H0 who has been in the conjoint study of 'never-subscribers'
ind_cj_never <- intersect(pop_H0$uido, data_nonsC_never$subjects$uido)

# Set corresponding indicators based on the above
pop_H0 <- pop_H0 %>% mutate(source = replace(source, uido %in% ind_sub, 1))
pop_H0 <- pop_H0 %>% mutate(source = replace(source, uido %in% ind_cj_never, 2))

# Find id's from data sources defined above

# Initialize u_ind
pop_H0$u_ind <- 0

# Case 1)
# Which id's from pop_H0 are in subscription data (i.e. earlier customers)
i1 <- pop_H0 %>%
  filter(source == 1) %>%
  select(uido) %>%
  pull()
# Search corresponding rows with the same original id from subscription data
i2 <- which(data_subs$subjects$uido %in% i1)
# Collect scaled id's from those rows
id_scaled <- data_subs$subjects[i2, "uid"] %>% pull()
# Save scaled id's above into the data set
pop_H0 <- pop_H0 %>% mutate(u_ind = replace(u_ind, uido %in% i1, id_scaled))


# Case 2)
# Which id's from pop_H0 are in 'never subscribed' conjoint data
i1 <- pop_H0 %>%
  filter(source == 2) %>%
  select(uido) %>%
  pull()
# Search corresponding rows with the same original id from conjoint data ('never subscribers')
i2 <- which(data_nonsC_never$subjects$uido %in% i1)
# Collect scaled id's from those rows
id_scaled <- data_nonsC_never$subjects[i2, "uid"] %>% pull()
# Save scaled id's above into the data set
pop_H0 <- pop_H0 %>% mutate(u_ind = replace(u_ind, uido %in% i1, id_scaled))

# Current customers

# The data source is always 1 (=earlier customers) because all individuals in pop_H1
# are included in the subscription data

pop_H1$source <- 1

# Find id's from data source defined above

# Initialize u_ind
pop_H1$u_ind <- 0

# All id's from pop_H1 are in subscription data
i1 <- pop_H1 %>%
  select(uido) %>%
  pull()

# Search corresponding rows with the same original id from subscription data
i2 <- which(data_subs$subjects$uido %in% i1)
# Collect scaled id's from those rows
id_scaled <- data_subs$subjects[i2, "uid"] %>% pull()
# Save scaled id's above into the data set
pop_H1 <- pop_H1 %>% mutate(u_ind = replace(u_ind, uido %in% i1, id_scaled))

# Add +1 for demographics for predictive simulation purposes
pop_H0 %<>% mutate(a = a + 1, g = g + 1, l = l + 1)
pop_H1 %<>% mutate(a = a + 1, g = g + 1, l = l + 1)


# Posterior predictive sampling -------------------------------------------

# Function for simulating the posterior predictive distribution
post_pred_sim <- function(posterior, pop_H0, pop_H1, price, N) {
  parameters <- colnames(posterior)
  N_A <- pop_H0$a |> n_distinct()
  N_G <- pop_H0$g |> n_distinct()
  N_L <- pop_H0$l |> n_distinct()

  # Posterior samples
  beta_const <- posterior[, "beta_const"]
  beta_age <- posterior[, grepl("beta_age", parameters)]
  beta_gender <- posterior[, grepl("beta_gender", parameters)]
  beta_location <- posterior[, grepl("beta_location", parameters)]
  alpha1 <- posterior[, "alpha1"]
  alpha2 <- numeric(N)
  alpha3 <- numeric(N)
  if ("alpha2" %in% parameters) {
    alpha2 <- posterior[, "alpha2"]
  }
  if ("alpha3" %in% parameters) {
    alpha3 <- posterior[, "alpha3"]
  }
  u_tau <- posterior[, "u_tau"]
  u_sub <- posterior[, grepl("u_sub", parameters)]
  if (ncol(beta_age) == 0L) {
    beta_age <- matrix(0.0, N, N_A)
    beta_gender <- matrix(0.0, N, N_G)
    beta_location <- matrix(0.0, N, N_L)
  }
  # If using all the data
  if (use_all_data) {
    kappa <- posterior[, "kappa"]
    u_nonC_never <- posterior[, grepl("u_nonC_never", parameters)]
  }

  # Observed data sets
  age0 <- pop_H0$a
  gender0 <- pop_H0$g
  location0 <- pop_H0$l
  prds0 <- pop_H0$prds
  source0 <- pop_H0$source
  u_ind0 <- pop_H0$u_ind

  age1 <- pop_H1$a
  gender1 <- pop_H1$g
  location1 <- pop_H1$l
  prds1 <- pop_H1$prds
  u_ind1 <- pop_H1$u_ind

  # Threshold price
  tp0 <- 0
  tp1 <- 0

  # Purchase probability
  p0 <- matrix(0, nrow = nrow(pop_H0), ncol = N)
  p1 <- matrix(0, nrow = nrow(pop_H1), ncol = N)

  # Iterate over posterior samples
  for (i in seq_len(N)) {
    # Calculate reference prices
    if (use_all_data) {
      tp0 <- exp(
        beta_const[i] +
        beta_age[i, age0] +
        beta_gender[i, gender0] +
        beta_location[i, location0] +
        ifelse(
          source0 == 1, 
          u_sub[i, u_ind0],
          ifelse(
            source0 == 2, 
            u_nonC_never[i, u_ind0],
            rnorm(nrow(pop_H0), mean =  0, sd = u_tau[i])
          )
        )
      )
    } else { # using only sales data
      tp0 <- exp(
        beta_const[i] +
        beta_age[i, age0] +
        beta_gender[i, gender0] +
        beta_location[i, location0] +
        ifelse(
          source0 == 1, 
          u_sub[i, u_ind0], 
          rnorm(nrow(pop_H0), mean =  0, sd = u_tau[i])
        )
      )
    }

    tp1 <- exp(
      beta_const[i] +
      beta_age[i, age1] +
      beta_gender[i, gender1] +
      beta_location[i, location1] +
      u_sub[i, u_ind1]
    )

    # Calculate purchase probabilities
    p0[, i] <- plogis((tp0 - price) * alpha1[i] + alpha2[i] * log1p(prds0) + alpha3[i] * (prds0 == 0))
    p1[, i] <- plogis((tp1 - price) * alpha1[i] + alpha2[i] * log1p(prds1) + alpha3[i] * (prds1 == 0))
  }

  # Return predictions
  list(p0 = p0, p1 = p1)
}


# Posterior predictive summaries ------------------------------------------


# Function for calling posterior predictive simulation and forming results
post_pred_summaries <- function(posterior, pop_H0, pop_H1, N0, N1, true_opt,
                                prices, M_post = NULL, varCost, fixCost) {
  # Initialize:

  # Data frame for results
  purchases <- as.matrix(data.frame(
    price = prices,
    EP = NA,
    EP_lwr95 = NA,
    EP_upp95 = NA,
    RMSE = NA
  ))

  # Vector of profits
  profit <- vector(mode = "list", length = nrow(purchases))
  
  # Merge the chains
  if (is.list(posterior)) {
    posterior <- do.call("rbind", posterior)
  }
  
  if (is.null(M_post)) {
    M_post <- nrow(posterior)
  }

  # Iterate over prices
  for (i in 1:nrow(purchases)) {
    # Simulate new purchases
    sim <- post_pred_sim(
      posterior = posterior,
      pop_H0 = pop_H0,
      pop_H1 = pop_H1,
      price = purchases[i, "price"],
      N = M_post
    )

    # Calculate expected gross profits and their CIs
    profit[[i]] <- N0 * colMeans(sim$p0) * (purchases[i, "price"] - varCost) + 
      N1 * colMeans(sim$p1) * (purchases[i, "price"] - varCost) - fixCost

    purchases[i, "EP"] <- mean(profit[[i]])
    purchases[i, "EP_lwr95"] <- quantile(profit[[i]], p = 0.025)
    purchases[i, "EP_upp95"] <- quantile(profit[[i]], p = 0.975)
    purchases[i, "RMSE"] <- sqrt(mean((profit[[i]] - true_opt[i, "EP"])^2))
  }

  # Collect and return the results
  profits <- data.frame(profit = unlist(profit))
  profits$price <- rep(prices, each = M_post)
  profits$iteration <- rep(seq_len(M_post), length(prices))

  list(purchases = purchases, profits = profits)
}


# Calculate expected profits ----------------------------------------------


prices <- seq(14, 18, by = 0.25)

# Set the price range to be 14 to 18
real_optimum %<>% subset(price >= 14 & price <= 18)
names(real_optimum)[2] <- "EP"
real_optimum$true_value <- 1

summary_all_data <- post_pred_summaries(
  posterior = samples,
  pop_H0 = pop_H0,
  pop_H1 = pop_H1,
  N0 = nrow(pop_H0),
  N1 = nrow(pop_H1),
  true_opt = real_optimum,
  prices = prices,
  varCost = 5,
  fixCost = 0
)


# Plots -------------------------------------------------------------------


# Add an indicator of estimates (0) vs. true values (1)
summary_all_data$purchases %<>% as.data.frame() %>% mutate(true_value = 0)

# Combine estimates and the true optimums into one data frame
optimums <- bind_rows(summary_all_data$purchases, real_optimum) %>%
  mutate(true_value = factor(true_value))

# Optimal price plot
opt_price_plot <- ggplot(data = optimums) +
  geom_point(
    aes(x = ifelse(true_value == 1, price + 0.04, price - 0.04), y = EP, shape = true_value)
  ) +
  scale_shape_manual(
    values = c(16, 1), 
    labels = c("Estimated expected gross profit", "True expected gross profit")
  ) +
  geom_linerange(
    aes(x = price - 0.04, ymin = EP_lwr95, ymax = EP_upp95),
    optimums |> filter(true_value == 0),
    inherit.aes = FALSE
  ) +
  scale_x_continuous(n.breaks = 8) +
  scale_y_continuous(
    limits = c(12000, 14500), 
    expand = expansion(mult = c(0.0022, 0.0022))
  ) +
  xlab("Price") +
  ylab("Expected gross profit") +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "inside",
    legend.justification = c(0, 0), 
    legend.position.inside = c(0.05, 0.055),
    legend.title = element_blank(), 
    legend.text = element_text(size = 10),
    legend.background = element_rect(fill = "white", color = "black")
  )

opt_price_plot

# Distribution of the optimal price
opt_dist <- summary_all_data$profits |>
  group_by(iteration) |>
  summarise(opt_price = price[which.max(profit)]) |>
  group_by(opt_price) |>
  count() |>
  ungroup() |>
  mutate(prop = n / sum(n))

# Plot the above
opt_price_distr <- ggplot(data = opt_dist, aes(x = factor(opt_price), y = prop)) +
  scale_y_continuous(
    breaks = seq(0, 0.5, by = 0.05), 
    limits = c(0, 0.40),
    expand = expansion(mult = c(0.0022, 0.0022))
  ) +
  geom_col(fill = "darkgray") +
  ylab("Probability of being the optimal price") +
  xlab("Price") +
  theme_bw(base_size = 12) +
  theme(panel.grid.major = element_blank())

opt_price_distr

# Create combined plots of both scenarios here
combined_plots <- opt_price_plot + opt_price_distr

if (use_all_data) {
  ggsave(
    filename = paste0("optimal_price_combined", suffix, ".pdf"), 
    plot = combined_plots, height = 5.5, width = 12, scale = 0.8
  )
} else { # using only sales data
  ggsave(
    filename = "optimal_price_combined_only_sales.pdf", 
    plot = combined_plots, height = 5.5, width = 12, scale = 0.8
  )
}


# Summaries and diagnostics -----------------------------------------------


posterior_draws <- posterior::bind_draws(lapply(samples, posterior::as_draws), along = "chain")

# Summaries and diagnostics for non-individual-specific parameters
sub <- posterior::subset_draws(
  posterior_draws, 
  variable = c(colnames(samples[[1]])[1:13], "u_tau")
)
sub_summ <- posterior::summarize_draws(sub)
sub_summ

# Diagnostics for all parameters
conv_all <- posterior::summarize_draws(posterior_draws, posterior::default_convergence_measures())
max(c(sub_summ$rhat, conv_all$rhat), na.rm = TRUE)
min(c(sub_summ$ess_bulk, conv_all$ess_bulk), na.rm = TRUE)
min(c(sub_summ$ess_bulk, conv_all$ess_tail), na.rm = TRUE)
