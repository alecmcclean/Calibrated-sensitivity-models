####################################################
### Darfur bounded odds ratio with M
### Use "no M" estimates for point estimates
### This script just quantifies uncertainty using
### the bootstrap.
####################################################

set.seed(20250624)

### Load data and functions
load("../intermediate/darfur/darfur_data.RData")
source("0_bdd_odds_functions.R")

#####################################################
### Create table of ranges for all covariates
### which will be used when estimating M

cov_ranges_dat <- 
  data.frame(var = colnames(darfur_covariates),
             max = apply(darfur_covariates, 2, max),
             min = apply(darfur_covariates, 2, min)) 


#######################################
### Bootstrap estimates

Gammas <- seq(1, 2, length.out = 5) 

NUM_BOOTS <- 100
BOOT_SIZE <- 1000
boot_est <- expand.grid(Gamma = Gammas,
                        boot = 1:NUM_BOOTS)
boot_est$max_confounding <- NA; boot_est$max_COV <- NA
boot_est$upper_est <- NA; boot_est$lower_est <- NA

start_time <- Sys.time()
results <- vector("list", NUM_BOOTS)
for (BOOT in seq_len(NUM_BOOTS)) {

  cat("\nBoot: ", BOOT, "out of", NUM_BOOTS)
    
  #####################################
  ### Resample for this bootstrap
  
  ROWS <- sample(1:nrow(darfur), BOOT_SIZE, replace = T)
  
  boot_darfur <- darfur[ROWS,]
  treatment <- boot_darfur$directlyharmed
  outcome <- boot_darfur$peacefactor
  covariates <- darfur_covariates[ROWS,]
  
  ##################################
  ### Estimate bounds
  
  ### Dataset for storing IF value estimates
  if_ests <- data.frame()
  
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
    
    cat("\nBoot: ", BOOT, "is at fold ", FOLD)
    
    ### Create training and estimation data
    train_treatment <- treatment[folds != FOLD]
    train_outcome <- outcome[folds != FOLD]
    train_covariates <- covariates[folds != FOLD,]
    
    est_treatment <- treatment[folds == FOLD]
    est_outcome <- outcome[folds == FOLD]
    est_covariates <- covariates[folds == FOLD,]
    
    # In training data, estimate pihat logistic regression with no interactions
    prop_mod <-  glm(directlyharmed ~ -1 + ., 
                     dplyr::bind_cols("directlyharmed" = train_treatment, train_covariates),
                     family = binomial(link = "logit"))
    
    # Calculate Mhat
    max_COV <- NA
    max_confounding <- 0
    for (COV in colnames(train_covariates)) {
      COV_confounding <- NA
      
      if (!grepl("village", COV)) {
        # Safe coefficient lookup
        COV_coef <- coef(prop_mod)[COV]
        
        if (!is.na(COV_coef)) {
          range_diff <- cov_ranges_dat$max[cov_ranges_dat$var == COV] -
            cov_ranges_dat$min[cov_ranges_dat$var == COV]
          COV_confounding <- abs(COV_coef * range_diff)
        }
        
      } else if (COV == "village_abu_gamra") {
        # Handle factor contrasts (village indicators)
        village_coefs <- coef(prop_mod)[grepl("village", names(coef(prop_mod)))]
        
        if (length(village_coefs) > 0) {
          for (i in village_coefs) {
            for (j in village_coefs) {
              COV_confounding <- abs(i - j)
              if (!is.na(COV_confounding) && COV_confounding > max_confounding) {
                max_COV <- "village"
                max_confounding <- COV_confounding
              }
            }
          }
        }
        next
      }
      
      # Global update if applicable
      if (!is.na(COV_confounding) && COV_confounding > max_confounding) {
        max_COV <- COV
        max_confounding <- COV_confounding
      }
    }
    
    
    # For each Gamma, nonparametrically estimate nuisance functions using 
    # max_confounding as input
    for (gamma in Gammas) {
      
      ADJ_GAMMA <- exp(log(gamma) / max_confounding)
      
      ### Estimate thetas
      coef_theta1_low <- optim(par = rep(0, ncol(train_covariates)),
                               fn = theta_optim,
                               method = "BFGS",
                               y = train_outcome[train_treatment == 1],
                               X = train_covariates[train_treatment == 1, ],
                               Gamma = ADJ_GAMMA)$par
      
      coef_theta1_high <- optim(par = rep(0, ncol(train_covariates)),
                                fn = theta_optim,
                                method = "BFGS",
                                y = train_outcome[train_treatment == 1],
                                X = train_covariates[train_treatment == 1, ],
                                Gamma = 1 / ADJ_GAMMA)$par
      
      coef_theta0_low <- optim(par = rep(1, ncol(train_covariates)),
                               fn = theta_optim,
                               method = "BFGS",
                               y = outcome[train_treatment == 0],
                               X = covariates[train_treatment == 0, ],
                               Gamma = ADJ_GAMMA)$par
      
      coef_theta0_high <- optim(par = rep(1, ncol(train_covariates)),
                                fn = theta_optim,
                                method = "BFGS",
                                y = outcome[train_treatment == 0],
                                X = covariates[train_treatment == 0, ],
                                Gamma = 1 / ADJ_GAMMA)$par
      
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
                            Gamma = ADJ_GAMMA)$par
      
      coef_nu1_high <- optim(par = rep(1, ncol(train_covariates)),
                             fn = nu_optim,
                             method = "BFGS",
                             y = train_outcome[train_treatment == 1],
                             X = train_covariates[train_treatment == 1, ],
                             theta = theta1_high[train_treatment == 1],
                             Gamma = ADJ_GAMMA)$par
      
      coef_nu0_low <- optim(par = rep(1, ncol(train_covariates)),
                            fn = nu_optim,
                            method = "BFGS",
                            y = train_outcome[train_treatment == 0],
                            X = train_covariates[train_treatment == 0, ],
                            theta = theta1_low[train_treatment == 0],
                            Gamma = ADJ_GAMMA)$par
      
      coef_nu0_high <- optim(par = rep(1, ncol(train_covariates)),
                             fn = nu_optim,
                             method = "BFGS",
                             y = train_outcome[train_treatment == 0],
                             X = train_covariates[train_treatment == 0, ],
                             theta = theta1_high[train_treatment == 0],
                             Gamma = ADJ_GAMMA)$par
      
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
      
      if_ests <- if_ests%>% bind_rows(
        data.frame(if_low = if_low, if_high = if_high, Gamma = gamma)
      )
    }
  }
  
  # Calculate lower bound and lower CI
  for (gamma in Gammas) {
    
    temp <- if_ests %>% filter(Gamma == gamma)
    
    # Calculate lower bound and lower CI
    boot_est$lower_est[boot_est$Gamma == gamma & boot_est$boot == BOOT] <- mean(temp$if_low)
    boot_est$lower_ci_lb[boot_est$Gamma == gamma & boot_est$boot == BOOT] <- 
      mean(temp$if_low) - 1.96 * sqrt(var(temp$if_low) / nrow(temp))
    
    # Calculate upper bound and upper CI
    boot_est$upper_est[boot_est$Gamma == gamma & boot_est$boot == BOOT] <- mean(temp$if_high)
    boot_est$upper_ci_ub[boot_est$Gamma == gamma & boot_est$boot == BOOT] <- 
      mean(temp$if_high) + 1.96 * sqrt(var(temp$if_high) / nrow(temp))
    
  }
  
  # Calculate max confounding on whole dataset
  prop_mod <-  glm(directlyharmed ~ -1 + ., 
                   dplyr::bind_cols("directlyharmed" = treatment, covariates),
                   family = binomial(link = "logit"))
  
  boot_est$max_confounding[boot_est$boot == BOOT] <- 0
  
  for (COV in colnames(covariates)) {
    COV_confounding <- NA
    
    if (!grepl("village", COV)) {
      COV_coef <- coef(prop_mod)[COV]
      
      if (!is.na(COV_coef)) {
        range_diff <- cov_ranges_dat$max[cov_ranges_dat$var == COV] -
          cov_ranges_dat$min[cov_ranges_dat$var == COV]
        COV_confounding <- abs(COV_coef * range_diff)
      }
      
    } else if (COV == "village_abu_gamra") {
      village_coefs <- coef(prop_mod)[grepl("village", names(coef(prop_mod)))]
      
      if (length(village_coefs) > 0) {
        for (i in village_coefs) {
          for (j in village_coefs) {
            COV_confounding <- abs(i - j)
            if (!is.na(COV_confounding) &&
                COV_confounding > boot_est$max_confounding[boot_est$boot == BOOT][1]) {
              boot_est$max_COV[boot_est$boot == BOOT] <- "village"
              boot_est$max_confounding[boot_est$boot == BOOT] <- COV_confounding
            }
          }
        }
      }
      next
    }
    
    if (!is.na(COV_confounding) &&
        COV_confounding > boot_est$max_confounding[boot_est$boot == BOOT][1]) {
      boot_est$max_COV[boot_est$boot == BOOT] <- COV
      boot_est$max_confounding[boot_est$boot == BOOT] <- COV_confounding
    }
  }
  
  results[[BOOT]] <- boot_est
}

# End timing
end_time <- Sys.time()

# Calculate elapsed time
elapsed_time <- end_time - start_time
print(elapsed_time)

### Pull results from list
boot_est <- data.frame()
for (i in 1:length(results)) {
  boot_est <- boot_est %>% bind_rows(results[[i]] %>% filter(!is.na(max_confounding), boot == i))
}


#####################
### Save

saveRDS(boot_est, file = "../intermediate/darfur/bdd_odds_bootstraps.RDS")

rm(list = ls(all = T))
gc()

