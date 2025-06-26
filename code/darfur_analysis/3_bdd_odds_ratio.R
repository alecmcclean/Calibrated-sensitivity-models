####################################################
### Darfur bounded odds ratio without M
####################################################

set.seed(20250624)

load("../intermediate/darfur/darfur_data.RData")
source("0_bdd_odds_functions.R")

#####################################################
### Create table of ranges for all covariates
### which will be used when estimating M

cov_ranges_dat <- 
  data.frame(var = colnames(darfur_covariates),
             max = apply(darfur_covariates, 2, max),
             min = apply(darfur_covariates, 2, min)) 


##########################################################
### Point estimate and CI without uncertainty in M

Gammas <- seq(1, 2, length.out = 5) 

### Dataset for final results
estimates <- expand.grid(Gamma = Gammas)
estimates$max_confounding <- 0; estimates$max_COV <- NA
estimates$upper_est <- NA; estimates$lower_est <- NA
estimates$upper_ci_ub <- NA; estimates$lower_ci_lb <- NA

### Dataset for storing IF value estimates
if_ests <- data.frame()

treatment <- darfur$directlyharmed
outcome <- darfur$peacefactor
covariates <- darfur_covariates

### Sample split
# Define the number of folds
NUM_FOLDS <- 5

# Shuffle the indices of the data
count_per_number <- nrow(covariates) %/% NUM_FOLDS

# Create a vector with the desired counts of each number
vector <- rep(1:NUM_FOLDS, each = count_per_number)

# If there are any remaining observations, distribute them randomly among the numbers
remaining <- nrow(covariates) %% NUM_FOLDS
if (remaining > 0) {
  vector <- c(vector, sample(1:NUM_FOLDS, remaining, replace = TRUE))
}

# Scramble to get the folds
folds <- sample(vector)

for (FOLD in 1:NUM_FOLDS) {
  
  ### Create training and estimation data
  train_treatment <- treatment[folds != FOLD]
  train_outcome <- outcome[folds != FOLD]
  train_covariates <- covariates[folds != FOLD,]
  
  est_treatment <- treatment[folds == FOLD]
  est_outcome <- outcome[folds == FOLD]
  est_covariates <- covariates[folds == FOLD,]
  
  # For each Gamma, nonparametrically estimate nuisance functions 
  for (gamma in Gammas) {
    
    cat("\nGamma: ", gamma)
    
    ### Estimate thetas
    coef_theta1_low <- optim(par = rep(1, ncol(train_covariates)),
                             fn = theta_optim,
                             method = "BFGS",
                             y = train_outcome[train_treatment == 1],
                             X = train_covariates[train_treatment == 1, ],
                             Gamma = gamma)$par
    
    coef_theta1_high <- optim(par = rep(1, ncol(train_covariates)),
                              fn = theta_optim,
                              method = "BFGS",
                              y = train_outcome[train_treatment == 1],
                              X = train_covariates[train_treatment == 1, ],
                              Gamma = 1 / gamma)$par
    
    coef_theta0_low <- optim(par = rep(1, ncol(train_covariates)),
                             fn = theta_optim,
                             method = "BFGS",
                             y = outcome[train_treatment == 0],
                             X = covariates[train_treatment == 0, ],
                             Gamma = gamma)$par
    
    coef_theta0_high <- optim(par = rep(1, ncol(train_covariates)),
                              fn = theta_optim,
                              method = "BFGS",
                              y = outcome[train_treatment == 0],
                              X = covariates[train_treatment == 0, ],
                              Gamma = 1 / gamma)$par
    
    ### Create theta vectors on training sample to estimate nu
    theta1_low <- as.vector(as.matrix(train_covariates) %*% as.matrix(coef_theta1_low))
    theta1_high <- as.vector(as.matrix(train_covariates) %*% as.matrix(coef_theta1_high))
    theta0_low <- as.vector(as.matrix(train_covariates) %*% as.matrix(coef_theta0_low))
    theta0_high <- as.vector(as.matrix(train_covariates) %*% as.matrix(coef_theta0_high))
    
    
    ### Estimate nus
    coef_nu1_low <- optim(par = rep(1, ncol(train_covariates)),
                          fn = nu_optim,
                          method = "BFGS",
                          y = train_outcome[train_treatment == 1],
                          X = train_covariates[train_treatment == 1, ],
                          theta = theta1_low[train_treatment == 1],
                          Gamma = gamma)$par
    
    coef_nu1_high <- optim(par = rep(1, ncol(train_covariates)),
                           fn = nu_optim,
                           method = "BFGS",
                           y = train_outcome[train_treatment == 1],
                           X = train_covariates[train_treatment == 1, ],
                           theta = theta1_high[train_treatment == 1],
                           Gamma = gamma)$par
    
    coef_nu0_low <- optim(par = rep(1, ncol(train_covariates)),
                          fn = nu_optim,
                          method = "BFGS",
                          y = train_outcome[train_treatment == 0],
                          X = train_covariates[train_treatment == 0, ],
                          theta = theta1_low[train_treatment == 0],
                          Gamma = gamma)$par
    
    coef_nu0_high <- optim(par = rep(1, ncol(train_covariates)),
                           fn = nu_optim,
                           method = "BFGS",
                           y = train_outcome[train_treatment == 0],
                           X = train_covariates[train_treatment == 0, ],
                           theta = theta1_high[train_treatment == 0],
                           Gamma = gamma)$par
    
    ### Estimate nu in estimation sample
    nu1_low <- as.vector(as.matrix(est_covariates) %*% as.matrix(coef_nu1_low))
    nu1_high <- as.vector(as.matrix(est_covariates) %*% as.matrix(coef_nu1_high))
    nu0_low <- as.vector(as.matrix(est_covariates) %*% as.matrix(coef_nu0_low))
    nu0_high <- as.vector(as.matrix(est_covariates) %*% as.matrix(coef_nu0_high)) 
    
    ### Estimate thetas in estimation sample
    theta1_low <- as.vector(as.matrix(est_covariates) %*% as.matrix(coef_theta1_low))
    theta1_high <- as.vector(as.matrix(est_covariates) %*% as.matrix(coef_theta1_high))
    theta0_low <- as.vector(as.matrix(est_covariates) %*% as.matrix(coef_theta0_low))
    theta0_high <- as.vector(as.matrix(est_covariates) %*% as.matrix(coef_theta0_high))
    
    
    ### Estimate propensity score nonparametrically
    prop_mod_np <- ranger::ranger(
      directlyharmed ~ -1 + .,
      dplyr::bind_cols("directlyharmed" = train_treatment, train_covariates))
    
    pihat_np <- predict(prop_mod_np, 
                        data = bind_cols("directly_harmed" = est_treatment,
                                         est_covariates))$predictions
    
    ### Calculate influence function values
    if_low <- as.vector(
      CalculateIFlow(est_treatment, est_outcome, theta1_low, theta0_high,
                     nu1_low, nu0_high, pihat_np, gamma)
    )
    
    if_high <- as.vector(
      CalculateIFhigh(est_treatment, est_outcome, theta1_high, theta0_low,
                      nu1_high, nu0_low, pihat_np, gamma)
    )
    
    if_ests <- if_ests %>% bind_rows(
      data.frame(if_low = if_low, if_high = if_high, Gamma = gamma)
    )
  }
}

# Calculate lower bound and lower CI
for (gamma in Gammas) {
  
  temp <- if_ests %>% filter(Gamma == gamma)
  
  # Calculate lower bound and lower CI
  estimates$lower_est[estimates$Gamma == gamma] <- mean(temp$if_low)
  estimates$lower_ci_lb[estimates$Gamma == gamma] <- 
    mean(temp$if_low) - 1.96 * sqrt(var(temp$if_low) / nrow(temp))
  
  # Calculate upper bound and upper CI
  estimates$upper_est[estimates$Gamma == gamma] <- mean(temp$if_high)
  estimates$upper_ci_ub[estimates$Gamma == gamma] <- 
    mean(temp$if_high) + 1.96 * sqrt(var(temp$if_high) / nrow(temp))
  
}

# Calculate max confounding on whole dataset
prop_mod <-  glm(directlyharmed ~ -1 + ., 
                 dplyr::bind_cols("directlyharmed" = treatment, covariates),
                 family = binomial(link = "logit"))

for (COV in colnames(covariates)) {
  if (!grepl("village", COV)) {
    COV_coef <- coef(prop_mod)[names(coef(prop_mod)) == COV]
    COV_confounding <- 
      abs(COV_coef * (cov_ranges_dat$max[cov_ranges_dat$var == COV] -
                        cov_ranges_dat$min[cov_ranges_dat$var == COV]))
    
  } else if (COV == "village_abu_gamra") {
    # check_max <- 0
    COV_coef <- coef(prop_mod)[grepl("village", names(coef(prop_mod)))]
    for (i in COV_coef) {
      for (j in COV_coef) {
        COV_confounding <- abs(i - j)
        if(COV_confounding > estimates$max_confounding[1]) {
          estimates$max_COV <- "village"
          estimates$max_confounding <- COV_confounding
        }
      }
    }
    
  } else {
    break
  }
  
  if(COV_confounding > estimates$max_confounding[1]) {
    estimates$max_COV <- COV
    estimates$max_confounding <- COV_confounding
  }
}

saveRDS(estimates, "../intermediate/darfur/bdd_odds_no_M.RDS")
rm(list = ls(all = T))
gc()
