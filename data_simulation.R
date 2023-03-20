
library(R6causal)
library(data.table)
library(dplyr)
library(magrittr)

rm(list = ls())
set.seed(2492021)


########################################################################################
#   GET DEMOGRAPHICS DATA AND CREATE POPULATION
########################################################################################


# Joint distribution of age, gender and location.
# "Vaesto asuinpaikan kaupunki-maaseutu-luokituksen seka sukupuolen ja ian mukaan, 2000-2020 , Tilastokeskus.
# Aineisto on ladattu Tilastokeskuksen rajapintapalvelusta 27.9.2021 lisenssilla CC BY 4.0"

finnAGL <- read.csv("finns_agl.csv", header = TRUE, sep = ";")

# Subset & classification of ages
finnAGL <- bind_cols(Ika = c(1, 2, 3, 4), bind_rows(
  finnAGL %>%
    slice(1:13) %>%
    select(-1) %>%
    colSums(),
  finnAGL %>%
    slice(14:28) %>%
    select(-1) %>%
    colSums(),
  finnAGL %>%
    slice(29:43) %>%
    select(-1) %>%
    colSums(),
  finnAGL %>%
    slice(44:58) %>%
    select(-1) %>%
    colSums())) %>% as.data.frame()

# Variables in the demographics data:
#
# ages:
# 0: 18-30
# 1: 31-45
# 2: 46-60
# 3: 61-75
#
# gender:
# 0: male
# 1: female
#
# location:
# 0: urban
# 1: rural
#

# Frequency table of the demographics combinations:
agl <- data.frame(expand.grid(
  agegroups = 1:4,
  gender = c("male", "female"),
  location = c("city", "countryside"),
  count = NA
))

# Calculate counts for each row of the table above
agl$count <- c(finnAGL$MiehetK, finnAGL$NaisetK, finnAGL$MiehetMS, finnAGL$NaisetMS)


############################################
#   Create population
############################################

# Create population dataframe:
#   * g = gender
#   * a = age
#   * l = location
population <- data.frame(id = integer(),
                         a = integer(),
                         g = factor(),
                         l = factor(),
                         prds = integer(),
                         urp = numeric(),
                         uk = numeric())

# Create every single unit in the population
for (i in 1:nrow(agl)) {
  population %<>% add_row(a = rep(agl[i,"agegroups"], agl[i,"count"]),
                          g = rep(agl[i,"gender"], agl[i,"count"]),
                          l = rep(agl[i,"location"], agl[i,"count"]))
}


# Define the size of population:
N <- nrow(population)

# Add following to the population dataframe:
#   * id = customer id (ido = will not be scaled)
#   * prds = lenght of recent subscription period
#   * urp = noise for reference prices of consumers
#   * uk = unobserved conjoint effect for consumers
population$id <- 1:N
population$prds <- 0
population$urp <- rnorm(N, 0, 0.1)
population$uk <- 0.75

# for demographics: substract 1 to get 0-level for purchase simulations
population$g <- as.numeric(population$g) - 1
population$a <- population$a - 1
population$l <- as.numeric(population$l) - 1


############################################
#   Define coefficients
############################################

# Define the constant and the effects of age, gender, and location
# on reference prices
b_0 <- 2.8
b_age1 <- -0.015
b_gen1 <- 0.01
b_loc1 <- -0.02

# Define the effects of reference price and periods on purchase probability
b_rpx <- 0.35
b_prds <- 0.45
b_first <- -0.3

# Set price increase
p_inc <- 0.5


##########################################################################################
#   SUBSCRIPTION SIMULATIONS
##########################################################################################

subscriptions <- SCM$new("subscriptions",
                         
                         uflist = list(
                           
                           # Purchase time
                           utime = "n : i", # from global variable
                           
                           # Customer id
                           uid = "n : id", # from global variable
                           
                           # Customer id (original i.e. will not be scaled)
                           uido = "n : id", # from global variable
                           
                           # Unobserved noise in reference price
                           urp = function(n) {
                             tp <- vector(length = n)
                             for (i in 1:n) {
                               tp[i] <- population[id[i], "urp"]
                             }
                             tp
                           }
                         ),
                         
                         vflist = list(
                           
                           # Get subscription period from 'population'
                           prds = function(uid) {
                             n_prds <- vector(length = length(uid))
                             for (i in 1:length(uid)) {
                               n_prds[i] <- population[uid[i], "prds"]
                             }
                             n_prds
                           },
                           
                           # Get age from 'population'
                           a = function(uid) {
                             age <- vector(length = length(uid))
                             for (i in 1:length(uid)) {
                               age[i] <- population[uid[i], "a"]
                             }
                             age
                           },
                           
                           # Get gender from 'population'
                           g = function(uid) {
                             gen <- vector(length = length(uid))
                             for (i in 1:length(uid)) {
                               gen[i] <- population[uid[i], "g"]
                             }
                             gen
                           },
                           
                           # Get location from 'population'
                           l = function(uid) {
                             loc <- vector(length = length(uid))
                             for (i in 1:length(uid)) {
                               loc[i] <- population[uid[i], "l"]
                             }
                             loc
                           },
                           
                           # Set price for the product concept:
                           x1 = function(utime) {
                             price1 <- vector(length = length(utime))
                             price1 <- p1
                             
                             if (utime[i] > 6) {
                               price1 <- (p1 + p_inc)
                               price1
                             }
                             if (utime[i] > 18) {
                               price1 <- (p1 + p_inc*2)
                               price1
                             }
                             if (utime[i] > 24) {
                               price1 <- X_prices
                             }
                             price1
                           },
                           
                           # Reference price
                           tp1 = "a, g, l, urp : exp(b_0 + b_age1*a + b_gen1*g + b_loc1*l + urp)",
                           
                           # Purchase probability
                           p = function(uid, x1, tp1, prds) {
                             purch_prob <- exp(b_rpx*(tp1 - x1) + b_prds*log(prds+1) + b_first*(prds==0))/
                               (1+exp(b_rpx*(tp1 - x1) + b_prds*log(prds+1) + b_first*(prds==0)))
                             purch_prob
                           },
                           
                           # Purchase choice
                           y = function(uid, p) {
                             decision <- rbinom(length(uid), 1, p)
                             decision
                           }
                         )
)



############################################
#   First round simulation (month = 1):
############################################

# Population size is N (initialized above)

# Total months to be simulated
M <- 24

# Number of consumers into purchase choice
n <- 1000

# Price of the product
p1 <- 16

# Simulation iteration round
i <- 1

# Sampling consumers into first purchase
id <- sample(population$id, size = n, replace = FALSE)

# Simulate purchases at time t=1
subscriptions$simulate(n)

# Current simulated data
dataCurrent <- subscriptions$simdata %>% as.data.frame()

# Update subscription periods
population[dataCurrent[dataCurrent$y == 1, "uid"], "prds"] <- population[dataCurrent[dataCurrent$y == 1, "uid"], "prds"] + 1
population[dataCurrent[dataCurrent$y == 0, "uid"], "prds"] <- 0 # reset if no purchase/cancellation

# Combine data into a larger data set (total sales)
data_subs <- dataCurrent

# Number of subscribed consumers at period t=1
n_cust <- sum(dataCurrent$y != 0)


############################################
#   Following simulations (periods t=1,2,...M)
############################################

for (i in 2:M) {
  
  # Update population for sampling
  n_noncust <- N - n_cust
  
  # New share of population into simulation
  m <- 1000
  
  # Update n
  n <- n_cust + m
  
  # Those who purchased at previous time
  customers_id <- dataCurrent[dataCurrent$y == 1, "uid"]
  
  # Sample id's for new potential customers
  non_customers_id <- sample(setdiff(population$id, customers_id), size = m, replace = FALSE)
  
  # Collect id's
  id <- c(customers_id, non_customers_id)
  
  # Simulate i_th subscriptions
  subscriptions$simulate(n)
  
  # Current simulation data
  dataCurrent <- subscriptions$simdata %>% as.data.frame()
  
  # Update the number of subscription periods
  population[dataCurrent[dataCurrent$y == 1, "uid"], "prds"] <- population[dataCurrent[dataCurrent$y == 1, "uid"], "prds"] + 1
  population[dataCurrent[dataCurrent$y == 0, "uid"], "prds"] <- 0 # reset if no purchase/cancellation
  
  # Combine datas
  data_subs <- unique(bind_rows(data_subs, dataCurrent))
  
  # Number of customers at t=i period (i.e. who subscribed)
  n_cust <- sum(dataCurrent$y != 0)
}

# Arrange data by id
data_subs %<>% arrange(uid)


##########################################################################################
#   CONJOINT STUDY SIMULATIONS
##########################################################################################

# Number of individuals to be picked up into a conjoint
n_units <- 200

# Number of tasks shown per customer
kC <- 10

# Set prices for the conjoint study
current_price <- max(data_subs$x1)
pricesC <- c(seq(current_price - 5, current_price + 5, by = 0.5))


# Sample a set for the conjoint of current customers:
id_current_cust <- data_subs %>% filter(utime == 24 & y == 1) %>% select(uid) %>% pull()
cj_subpop_subs <- population %>% filter((id %in% id_current_cust)) %>% sample_n(size = n_units, replace = FALSE)


# Sample a set for the conjoint of earlier customers:
id_non_cust2 <- setdiff(data_subs$uid, id_current_cust)
id_non_cust <- vector() # save id's here

# collect those id's with at least one subscription in the purchase history
for(i in 1:length(id_non_cust2)) {
  if(data_subs %>% filter(uid == id_non_cust2[i]) %>% select(y) %>% sum != 0) {
    id_non_cust <- c(id_non_cust, id_non_cust2[i])
  }
}

cj_subpop_nons <- population %>% filter((id %in% id_non_cust)) %>% sample_n(size = n_units, replace = FALSE)


# Sample a set for the conjoint of never been customers:
id_non_cust_never <- population %>% filter(!id %in% union(id_current_cust, id_non_cust)) %>% select(id) %>% pull()
cj_subpop_nons_never <- population %>% filter((id %in% id_non_cust_never)) %>% sample_n(size = n_units, replace = FALSE)


#############################################
#   Simulations for Non-subscribers (never been a customer)
#############################################

# Create conjoint frames
cj_subpop_nons_never <- do.call("rbind", replicate(n = kC, cj_subpop_nons_never, simplify = FALSE))
# Order data_subs by id
cj_subpop_nons_never %<>% arrange(id) 
# Add a column including information of price combinations to be shown
cj_subpop_nons_never$ux <- as.vector(replicate(n = n_units, expr = sample(1:length(pricesC), kC, replace = FALSE)))

# Conjoint simulation
conjoint_sim_nons_never <- SCM$new("conjoint_sim_nons_never",
                                   
                                   uflist = list(
                                     
                                     # Number of conjoint query
                                     utime = function(n) {
                                       rep(1:10, n_units)
                                     },
                                     
                                     # Collect id
                                     uid = "n : cj_subpop_nons_never$id",
                                     
                                     # Collect id (original i.e. will not be scaled)
                                     uido = "n : cj_subpop_nons_never$id",
                                     
                                     # Collect price combination to be shown
                                     ux = "n : cj_subpop_nons_never$ux",
                                     
                                     # Unobserved noise in reference price
                                     urp = "n : cj_subpop_nons_never$urp"
                                   ),
                                   
                                   vflist = list(
                                     
                                     # Subscription periods so far
                                     prds = function(uid) {
                                       n_prds <- vector(length = length(uid))
                                       for (i in 1:length(uid)) {
                                         n_prds[i] <- population[uid[i], "prds"]
                                       }
                                       n_prds
                                     },
                                     
                                     # Get age
                                     a = "uid : cj_subpop_nons_never$a",
                                     
                                     # Get gender
                                     g = "uid : cj_subpop_nons_never$g",
                                     
                                     # Get location
                                     l = "uid : cj_subpop_nons_never$l",
                                     
                                     # Get conjoint effect
                                     k = "uid : cj_subpop_nons_never$uk",
                                     
                                     # Prices to be shown
                                     x1 = function(ux) {
                                       x <- vector(length = length(ux))
                                       for (i in 1:length(ux)) {
                                         x[i] <- pricesC[ux[i]]
                                       }
                                       x
                                     },
                                     
                                     # Reference price
                                     tp1 = "a, g, l, urp : exp(b_0 + b_age1*a + b_gen1*g + b_loc1*l + urp)",
                                     
                                     # Choice probability
                                     p = function(uid, x1, tp1, prds, k) {
                                       choice_prob <- exp(b_rpx*(tp1 + k - x1) + b_prds*log(prds+1) + b_first*(prds==0))/
                                         (1+exp(b_rpx*(tp1 + k - x1) + b_prds*log(prds+1) + b_first*(prds==0)))
                                       choice_prob
                                     },
                                     
                                     # Choice
                                     y = function(uid, p) {
                                       decision <- rbinom(length(uid), 1, p)
                                       decision
                                     }
                                   )
)

# Simulate (non-subscribers)
conjoint_sim_nons_never$simulate(n_units * kC)
data_nonsC_never <- conjoint_sim_nons_never$simdata %>% as.data.frame()


#############################################
#   Simulations for Non-subscribers (earlier customers)
#############################################

# Create conjoint frames
cj_subpop_nons <- do.call("rbind", replicate(n = kC, cj_subpop_nons, simplify = FALSE))
# Order data_subs by id
cj_subpop_nons %<>% arrange(id) 
# Add a column including information of price combinations to be shown
cj_subpop_nons$ux <- as.vector(replicate(n = n_units, expr = sample(1:length(pricesC), kC, replace = FALSE)))

# Conjoint simulation
conjoint_sim_nons <- SCM$new("conjoint_sim_nons",
                             
                             uflist = list(
                               
                               # Number of conjoint query
                               utime = function(n) {
                                 rep(1:10, n_units)
                               },
                               
                               # Collect id
                               uid = "n : cj_subpop_nons$id",
                               
                               # Collect id (original i.e. will not be scaled)
                               uido = "n : cj_subpop_nons$id",
                               
                               # Collect price combination to be shown
                               ux = "n : cj_subpop_nons$ux",
                               
                               # Unobserved noise in reference price
                               urp = "n : cj_subpop_nons$urp"
                             ),
                             
                             vflist = list(
                               
                               # Subscription periods so far
                               prds = function(uid) {
                                 n_prds <- vector(length = length(uid))
                                 for (i in 1:length(uid)) {
                                   n_prds[i] <- population[uid[i], "prds"]
                                 }
                                 n_prds
                               },
                               
                               # Get age
                               a = "uid : cj_subpop_nons$a",
                               
                               # Get gender
                               g = "uid : cj_subpop_nons$g",
                               
                               # Get location
                               l = "uid : cj_subpop_nons$l",
                               
                               # Get conjoint effect
                               k = "uid : cj_subpop_nons$uk",
                               
                               # Prices to be shown
                               x1 = function(ux) {
                                 x <- vector(length = length(ux))
                                 for (i in 1:length(ux)) {
                                   x[i] <- pricesC[ux[i]]
                                 }
                                 x
                               },
                               
                               # Reference price
                               tp1 = "a, g, l, urp : exp(b_0 + b_age1*a + b_gen1*g + b_loc1*l + urp)",
                               
                               # Choice probability
                               p = function(uid, x1, tp1, prds, k) {
                                 choice_prob <- exp(b_rpx*(tp1 + k - x1) + b_prds*log(prds+1) + b_first*(prds==0))/
                                   (1+exp(b_rpx*(tp1 + k - x1) + b_prds*log(prds+1) + b_first*(prds==0)))
                                 choice_prob
                               },
                               
                               # Choice
                               y = function(uid, p) {
                                 decision <- rbinom(length(uid), 1, p)
                                 decision
                               }
                             )
)

# Simulate (non-subscribers)
conjoint_sim_nons$simulate(n_units * kC)
data_nonsC <- conjoint_sim_nons$simdata %>% as.data.frame()


#############################################
#   Simulations for subscribers (current customers)
#############################################

# Create conjoint frames
cj_subpop_subs <- do.call("rbind", replicate(n = kC, cj_subpop_subs, simplify = FALSE))
# Order data by id
cj_subpop_subs %<>% arrange(id)
# Add a column including information of price combinations to be shown
cj_subpop_subs$ux <- as.vector(replicate(n = n_units, expr = sample(1:length(pricesC), kC, replace = FALSE)))


# Conjoint simulation
conjoint_sim_subs <- SCM$new("conjoint_sim_subs",
                                    
                             uflist = list(
                               
                               # Number of conjoint query
                               utime = function(n) {
                                 rep(1:10, n_units)
                               },
                               
                               # Collect id
                               uid = "n : cj_subpop_subs$id",
                               
                               # Collect id (original i.e. will not be scaled)
                               uido = "n : cj_subpop_subs$id",
                               
                               # Collect price combination to be shown
                               ux = "n : cj_subpop_subs$ux",
                               
                               # Unobserved noise in reference price
                               urp = "n : cj_subpop_subs$urp"
                               
                             ),
                             
                             vflist = list(
                               
                               # Subscription periods so far
                               prds = function(uid) {
                                 n_prds <- vector(length = length(uid))
                                 for (i in 1:length(uid)) {
                                   n_prds[i] <- population[uid[i], "prds"]
                                 }
                                 n_prds
                               },
                               
                               # Get age
                               a = "uid : cj_subpop_subs$a",
                               
                               # Get gender
                               g = "uid : cj_subpop_subs$g",
                               
                               # Get location
                               l = "uid : cj_subpop_subs$l",
                               
                               # Get conjoint effect
                               k = "uid : cj_subpop_subs$uk",
                               
                               # Price combinations to be shown
                               x1 = function(ux) {
                                 x <- vector(length = length(ux))
                                 for (i in 1:length(ux)) {
                                   x[i] <- pricesC[ux[i]]
                                 }
                                 x
                               },
                               
                               # Reference price
                               tp1 = "a, g, l, urp : exp(b_0 + b_age1*a + b_gen1*g + b_loc1*l + urp)",
                               
                               # Choice probability
                               p = function(uid, x1, tp1, prds, k) {
                                 choice_prob <- exp(b_rpx*(tp1 + k - x1) + b_prds*log(prds+1) + b_first*(prds==0))/
                                   (1+exp(b_rpx*(tp1 + k - x1) + b_prds*log(prds+1) + b_first*(prds==0)))
                                 choice_prob
                               },
                               
                               # Choice
                               y = function(uid, p) {
                                 decision <- rbinom(length(uid), 1, p)
                                 decision
                               }
                             )
)

# Simulate (subscribers)
conjoint_sim_subs$simulate(n_units * kC)
data_subsC <- conjoint_sim_subs$simdata %>% as.data.frame()


##########################################################################################
#   FUTURE DATA SETS (t=25)
##########################################################################################


# Update time for simulating the next period/month
i <- 25

# Update n
n <- n_cust + m

# Those who purchased at previous time
customers_id <- dataCurrent[dataCurrent$y == 1, "uid"]

# Sample id for new potential customers (the rest of the whole population)
non_customers_id <- sample(setdiff(population$id, customers_id), size = m, replace = FALSE)

# Collect above id's
id <- c(customers_id, non_customers_id)


#############################################
#   Data for current customers
#############################################

# Data frame for those who subscribed at t=24
pop_H1 <- dataCurrent %>% filter(y==1) %>% as.data.frame()
pop_H1 %<>% select(utime, uid, uido, a, g, l, prds)
pop_H1$x1 <- NA; pop_H1$y <- NA # initialize
pop_H1$utime <- i # update time
pop_H1$prds <- pop_H1$prds + 1 # update the number of periods


#############################################
#   Data for earlier customers & never been customers
#############################################

# Find customers that have subscribed sometimes
non_customers_id1 <- intersect(non_customers_id, data_subs$uido)

# id_non_cust:
#non_customers_id

# Find corresponding data
data_H0_1 <- data_subs %>% filter(uid %in% non_customers_id1) %>% group_by(uid) %>% slice_tail()

# Select variables needed in new purchase simulations
data_H0_1 %<>% select(uid, uido, a, g, l, prds) %>% as.data.frame()

# Find customers that have never subscribed
non_customers_id2 <- setdiff(non_customers_id,non_customers_id1)

# Collect demographics from new customers and add to data frame
data_H0_2 <- population %>% filter(( population$id %in% non_customers_id2)) %>% as.data.frame()

# Rename columns for combining
data_H0_2 <- rename(data_H0_2, c(uid = id, uido = id))

# Clear scaled ID (uid)
data_H0_2$uid <- NA

# Select variables needed in new purchase simulations
data_H0_2 %<>% select(uid, uido, a, g, l, prds)

# Combine datas:
data_H0 <- rbind(data_H0_1, data_H0_2)

# Create a data frame
pop_H0 <- data.frame(utime = i, # update time
                     uid = data_H0$uid,
                     uido = data_H0$uido,
                     a = data_H0$a,
                     g = data_H0$g,
                     l = data_H0$l,
                     x1 = NA, # initialize
                     y = NA, # initialize
                     prds = 0) # reset the number of prds


# Save both t=25 datas
save(pop_H0, pop_H1, N, file = "newdata.RData")


##########################################################################################
#   MODIFY DATAS FOR MODELING
##########################################################################################


# 1) Subscription data:

# Separate purchases and customers
data_subs <- list(
  subjects = as_tibble(data_subs %>% select(utime, uid, uido, a, g, l, prds)) %>% arrange(uid),
  tasks = as_tibble(data_subs %>% select(utime, uid, uido, a, g, l, prds, x1, y)) %>% arrange(uid)
)
# Create subjects data
data_subs$subjects <- data_subs$subjects %>% distinct(uid, .keep_all = TRUE)


# 2) Conjoint data of current subscribers:

# Separate purchases and customers
data_subsC <- list(
  subjects = as_tibble(data_subsC %>% select(utime, uid, uido, a, g, l, prds)) %>% arrange(uid),
  tasks = as_tibble(data_subsC %>% select(utime, uid, uido, a, g, l, prds, x1, y)) %>% arrange(uid)
)
# Create subjects data
data_subsC$subjects <- data_subsC$subjects %>% distinct(uid, .keep_all = TRUE)


# 3) Conjoint data of earlier subscribers:

# Separate purchases and customers
data_nonsC <- list(
  subjects = as_tibble(data_nonsC %>% select(utime, uid, uido, a, g, l, prds)) %>% arrange(uid),
  tasks = as_tibble(data_nonsC %>% select(utime, uid, uido, a, g, l, prds, x1, y)) %>% arrange(uid)
)
# Create subjects data
data_nonsC$subjects <- data_nonsC$subjects %>% distinct(uid, .keep_all = TRUE)


# 4) Conjoint data of never subscribers:

# Separate purchases and customers
data_nonsC_never <- list(
  subjects = as_tibble(data_nonsC_never %>% select(utime, uid, uido, a, g, l, prds)) %>% arrange(uid),
  tasks = as_tibble(data_nonsC_never %>% select(utime, uid, uido, a, g, l, prds, x1, y)) %>% arrange(uid)
)
# Create subjects data
data_nonsC_never$subjects <- data_nonsC_never$subjects %>% distinct(uid, .keep_all = TRUE)


#############################################
#   ID scaling:  
#   for 2) and 3) above: find the corresponding id's in the subscription data
#############################################

# 2)
# Collect id's which corresponds customers in both subscription data and current subscribers' conjoint data

id_collector_subs <- vector() # initialize vector of ids

# Loop over conjoint data of subscribers
for (i in 1:length(data_subsC$subjects$uid)) {
  ind <- data_subsC$subjects$uid[i]
  if(ind %in% data_subs$subjects$uid) {
    id_collector_subs[i] <- which(ind == data_subs$subjects$uid)
  }
}

id_collector_subs <- na.omit(id_collector_subs) %>% as.vector()
# This gives us the rows in subscription data (subjects dataset)
# whose corresponding id is in subscribers' conjoint data (subjects dataset)


# 3)
# Collect id's which corresponds customers in both subscription data and non-subscribers' (earlier) conjoint data

id_collector_nonsubs <- vector() # initialize vector of ids

# Loop over conjoint data of non-subscribers
for (i in 1:length(data_nonsC$subjects$uid)) {
  ind <- data_nonsC$subjects$uid[i]
  if(ind %in% data_subs$subjects$uid) {
    id_collector_nonsubs[i] <- which(ind == data_subs$subjects$uid)
  }
}

id_collector_nonsubs <- na.omit(id_collector_nonsubs) %>% as.vector()
# This gives us the rows in subscription data (subjects dataset)
# whose corresponding id is in non-subscribers' (earlier) conjoint data (subjects dataset)


# Re-scale id's
data_subs$subjects %<>% mutate(uid = dense_rank(uid))
data_subs$tasks %<>% mutate(uid = dense_rank(uid)) 

data_subsC$subjects %<>% mutate(uid = dense_rank(uid))
data_subsC$tasks %<>% mutate(uid = dense_rank(uid))

data_nonsC$subjects %<>% mutate(uid = dense_rank(uid))
data_nonsC$tasks %<>% mutate(uid = dense_rank(uid))

data_nonsC_never$subjects %<>% mutate(uid = dense_rank(uid))
data_nonsC_never$tasks %<>% mutate(uid = dense_rank(uid))


#############################################
#   Add 1 for demographic variables for NIMBLE estimation purposes
#############################################

data_subs$subjects %<>% mutate(a=a+1, g=g+1, l=l+1)
data_subs$tasks %<>% mutate(a=a+1, g=g+1, l=l+1)

data_subsC$subjects %<>% mutate(a=a+1, g=g+1, l=l+1)
data_subsC$tasks %<>% mutate(a=a+1, g=g+1, l=l+1)

data_nonsC$subjects %<>% mutate(a=a+1, g=g+1, l=l+1)
data_nonsC$tasks %<>% mutate(a=a+1, g=g+1, l=l+1)

data_nonsC_never$subjects %<>% mutate(a=a+1, g=g+1, l=l+1)
data_nonsC_never$tasks %<>% mutate(a=a+1, g=g+1, l=l+1)


##########################################################################################
#   SAVE DATAS TO BE USED IN MODELING
##########################################################################################

save(data_subs, data_subsC, data_nonsC, data_nonsC_never, id_collector_subs, id_collector_nonsubs, population, N, file = "cabcData.RData")


##########################################################################################
#   FUTURE DATA PURCHASE SIMULATIONS (t=25) 
##########################################################################################

# Update time for simulating the next period/month
i <- 25

# Update n
n <- n_cust + m

# Define gross profit formula parameters:

# Variable cost of the product
varCost <- 5

# price range
X_prices_list <- seq(14,20, by=0.25)

# Initialize gross profit matrix
gross_profit <- vector(length = length(X_prices_list))

# Size of H0 and H1 datas
N0 <- nrow(pop_H0)
N1 <- nrow(pop_H1)

for(j in 1:length(X_prices_list)) {
  
  # Set a price
  X_prices <- X_prices_list[j]
  
  # Simulate purchases
  subscriptions$simulate(n)
  data_t25 <- subscriptions$simdata %>% as.data.frame()
  
  # Calculate purchase probabilities
  p0 <- data_t25 %>% filter(uido %in% pop_H0$uido) %>% select(p) %>% unlist() %>% mean()
  p1 <- data_t25 %>% filter(uido %in% pop_H1$uido) %>% select(p) %>% unlist() %>% mean()
  
  # Calculate gross profit
  gross_profit[j] <- N0*p0*(X_prices_list[j] - varCost) + N1*p1*(X_prices_list[j] - varCost)
  
}

head(gross_profit)

# Plot real optimums
real_optimum <- data.frame(price = X_prices_list, grossprofit = gross_profit)
plot(x = real_optimum$price, y = real_optimum$grossprofit, lwd=2, type = "o", col = 2)

# Save the expected gross profit at t=25
save(real_optimum, file = "real_optimum.RData")

