
library(dplyr)
library(magrittr)
library(ggplot2)
library(patchwork)

rm(list=ls())
load("cabcData.RData")
load("newdata.RData")

# indicator: FALSE means we use only sales data
use_all_data <- TRUE

ifelse(use_all_data==TRUE, load("post_samples.RData"), load("post_samples_only_sales.RData"))
set.seed(2492021)


########################################################################################
#   PREPARATIONS FOR THE t=25 DATA 
########################################################################################


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
i1 <- pop_H0 %>% filter(source==1) %>% select(uido) %>% pull
# Search corresponding rows with the same original id from subscription data
i2 <- which(data_subs$subjects$uido %in% i1)
# Collect scaled id's from those rows
id_scaled <- data_subs$subjects[i2,"uid"] %>% pull
# Save scaled id's above into the data set
pop_H0 <- pop_H0 %>% mutate(u_ind = replace(u_ind, uido %in% i1, id_scaled))


# Case 2)
# Which id's from pop_H0 are in 'never subscribed' conjoint data
i1 <- pop_H0 %>% filter(source==2) %>% select(uido)%>%pull
# Search corresponding rows with the same original id from conjoint data ('never subscribers')
i2 <- which(data_nonsC_never$subjects$uido %in% i1)
# Collect scaled id's from those rows
id_scaled <- data_nonsC_never$subjects[i2,"uid"] %>% pull
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
i1 <- pop_H1 %>% select(uido) %>% pull
# Search corresponding rows with the same original id from subscription data
i2 <- which(data_subs$subjects$uido %in% i1)
# Collect scaled id's from those rows
id_scaled <- data_subs$subjects[i2,"uid"] %>% pull
# Save scaled id's above into the data set
pop_H1 <- pop_H1 %>% mutate(u_ind = replace(u_ind, uido %in% i1, id_scaled))


############################################
# Add +1 for demographics for predictive simulation purposes
############################################

pop_H0 %<>% mutate(a=a+1,g=g+1,l=l+1)
pop_H1 %<>% mutate(a=a+1,g=g+1,l=l+1)


########################################################################################
#   POSTERIOR PREDICTIVE SAMPLING
########################################################################################

# Function for simulating posterior predictive
q1_sim <- function(posterior, posterior_u, pop_H0, pop_H1, price, M_post) {
  
  # Total number of posterior samples
  N <- posterior$beta_age %>% nrow()
  
  # Posterior samples
  beta_const <- posterior$beta_const
  beta_age <- posterior$beta_age
  beta_gender <- posterior$beta_gender
  beta_location <- posterior$beta_location
  alpha1 <- posterior$alpha1
  alpha2 <- posterior$alpha2
  alpha3 <- posterior$alpha3
  u_tau <- posterior$u_tau
  u_sub <- posterior_u$u_sub
  
  # If using all the data
  if(use_all_data == TRUE) {
    kappa <- posterior$kappa
    u_nonC_never <- posterior_u_nonC_never$u_nonC_never
  }
  
  # Set a number of samples to be used from posterior samples
  N <- M_post
  
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
  for (i in 1:N) {
    
    # Calculate reference prices
    if(use_all_data == TRUE) {
      tp0 <- exp(beta_const[i] +
                   beta_age[i,age0] +
                   beta_gender[i,gender0] +
                   beta_location[i,location0] +
                   ifelse(source0 == 1, u_sub[i, u_ind0],
                          ifelse(source0 == 2, u_nonC_never[i, u_ind0],
                                 rnorm(nrow(pop_H0), mean =  0, sd = u_tau[i]))))
    }else{ #using only sales data
      tp0 <- exp(beta_const[i] +
                   beta_age[i,age0] +
                   beta_gender[i,gender0] +
                   beta_location[i,location0] +
                   ifelse(source0 == 1, u_sub[i, u_ind0], rnorm(nrow(pop_H0), mean =  0, sd = u_tau[i])))
    }
    
    tp1 <- exp(beta_const[i] +
                 beta_age[i,age1] +
                 beta_gender[i,gender1] +
                 beta_location[i,location1] +
                 u_sub[i, u_ind1])
    
    
    # Calculate purchase probabilities
    p0[,i] <- plogis((tp0 - price)*alpha1[i] + alpha2[i]*log(prds0+1) + alpha3[i]*(prds0==0))
    
    p1[,i] <- plogis((tp1 - price)*alpha1[i] + alpha2[i]*log(prds1+1) + alpha3[i]*(prds1==0))
    
  }
  
  # Return predictions
  list(p0=p0, p1=p1)
}


########################################################################################
#   POSTERIOR PREDICTIVE SUMMARIES
########################################################################################

# Function for calling posterior predictive simulation and forming results
postPred_summaries <- function(pop_H0, pop_H1, N0, N1, prices, M_post, varCost, fixCost) {
  
  # Initialize:
  
  # Data frame for results
  purchases <- as.matrix(data.frame(expand.grid(price = prices,
                                                EP = NA,
                                                EP_lwr95 = NA,
                                                EP_upp95 = NA)))
  
  # Vector of profits
  profit <- vector(mode = "list", length = nrow(purchases))
  
  # Iterate over prices
  for(i in 1:nrow(purchases)) {
    
    # Simulate new purchases
    sim <- q1_sim(posterior = posterior,
                  posterior_u = posterior_u,
                  pop_H0 = pop_H0,
                  pop_H1 = pop_H1,
                  price = purchases[i, "price"],
                  M_post = M_post)
    
    # Calculate expected gross profits and their CIs
    profit[[i]] <- N0*(colMeans(sim$p0))*(purchases[i, "price"] - varCost) + N1*(colMeans(sim$p1))*(purchases[i, "price"] - varCost) - fixCost
    
    purchases[i,"EP"] <- mean(profit[[i]])
    purchases[i,"EP_lwr95"] <- quantile(profit[[i]],p=0.025)
    purchases[i,"EP_upp95"] <- quantile(profit[[i]],p=0.975)
    
  }
  
  # Collect and return the results
  profits <- data.frame(profit = unlist(profit))
  profits$price <- rep(prices, each = M_post)
  profits$iteration <- rep(seq_len(M_post), length(prices))
  
  list(purchases=purchases, profits=profits)
}


########################################################################################
#   CALCULATE EXPECTED PROFITS
########################################################################################

prices <- seq(14,18, by=0.25)

summaryAllData <- postPred_summaries(pop_H0 = pop_H0,
                                     pop_H1 = pop_H1,
                                     N0 = nrow(pop_H0),
                                     N1 = nrow(pop_H1) ,
                                     prices = prices,
                                     M_post = length(posterior$beta_const), # number of posterior samples used
                                     varCost = 5,
                                     fixCost = 0)

summaryAllData


########################################################################################
#   PLOT THE RESULTS
########################################################################################

# Load estimated expected gross profits
#load("optimal_prices.RData")

# Add an indicator of estimates (0) vs. true values (1)
summaryAllData$purchases %<>% as.data.frame() %>% mutate(true_value = 0)

# Load true expected gross profits
load("real_optimum.RData")
# Set the price range to be 14 to 18
real_optimum %<>% subset(price >= 14 & price <= 18)
names(real_optimum)[2] <- "EP"
real_optimum$true_value <- 1

# Combine estimates and the true optimums into one data frame
optimums <- bind_rows(summaryAllData$purchases, real_optimum) %>% 
  mutate(true_value = factor(true_value))

# Optimal price plot
opt_price_plot <- ggplot(data = optimums) +
  geom_point(aes(x = ifelse(true_value == 1, price + 0.04, price - 0.04), y = EP, shape = true_value)) +
  scale_shape_manual(values = c(16,1), labels =c("Estimated expected gross profit","True expected gross profit")) +
  geom_linerange(aes(x = price - 0.04, ymin = EP_lwr95, ymax = EP_upp95)) +
  scale_x_continuous(n.breaks=8) +
  xlab("Price") + ylab("Expected gross profit") +
  theme_bw() +
  theme(legend.position = c(0.4, 0.25), legend.title = element_blank(), legend.text=element_text(size=10),
        legend.background = element_rect(fill = "white", color = "black"))
  
opt_price_plot


# Distribution of the optimal price
opt_dist <- summaryAllData$profits |>
  group_by(iteration) |>
  summarise(opt_price = price[which.max(profit)]) |>
  group_by(opt_price) |>
  count() |>
  ungroup() |>
  mutate(prop = n/sum(n))


# Plot the above
opt_price_distr <- ggplot(data = opt_dist,
            aes(x = factor(opt_price), y = prop)) +
  scale_y_continuous(breaks = seq(0,0.5, by=0.05)) +
  geom_col(fill="darkgray") +
  ylab("Probability of being the optimal price") +
  xlab("Price") +
  theme_bw()
opt_price_distr


# Create combined plots of both scenarios here
combined_plots <- opt_price_plot + opt_price_distr


########################################################################################
#   SAVE THE RESULTS AND PLOTS
########################################################################################

# Save the results and plots
if(use_all_data == TRUE) {
  save(summaryAllData, opt_dist, file = "optimal_prices.RData")
  ggsave(filename = "optimal_price_combined.pdf", plot = combined_plots, scale = 0.7, height = 6)
}else{ #using only sales data
  save(summaryAllData, opt_dist, file = "optimal_prices_only_sales.RData")
  ggsave(filename = "optimal_price_combined_only_sales.pdf", plot = combined_plots, scale = 0.7, height = 6)
}

