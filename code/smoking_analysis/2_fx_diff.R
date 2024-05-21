####################################################
### Smoking effect differences analysis
####################################################

###################################
### Load data

load("../../intermediate/smoking/smoking_data.RData")

### Estimate ATE with full data 
full_est <- 
  npcausal::ate(y = smoking$dbirwt,
                a = smoking$is_smoke,
                x = select(smoking, -dbirwt, -is_smoke),
                nsplits = 5,
                sl.lib = c("SL.mean", "SL.ranger"))

full_ifvals <- full_est$ifvals$a1 - full_est$ifvals$a0
full_ate_est <- full_est$res$est[3]
full_ate_lb <- full_est$res$ci.ll[3]
full_ate_ub <- full_est$res$ci.ul[3]

# Print results
cat("ATE: ", full_ate_est, 
    "\n95% CI: [", full_ate_lb, ", ", full_ate_ub, "]")



#####################################
### Maximum leave-one-out analysis

### Keep track of covariate corresponding to max confounding
max_confounding <- 0
max_covariate <- "none"

COVS <- c("dmar", "race", "foreignb", "alcohol", "disllbu", "tripre",
          "adequac", "nprevisu", "ddeadkids", "dbirmon", "dcntyfipb",
          "age", "educ", "nprevist", "disllb",  "dlivord")

for (COV in COVS) {
  
  # Data cleaning --- remove all reasonable groups of variables as if they
  # were one variable (e.g., race, age, education of parents)
  if (COV == "race") {# Remove all race 
    
    COV <- c("mwhite", "mblack", "mhispan", "fwhite", "fblack", "fhispan")
    
  } else if (COV %in% c("age", "educ")) {
    
    COV <- colnames(smoking) %>% .[grepl(COV, .)]
    
    # Remove all binary variables that are levels of the same factor variables together
  } else if (COV %in% c("tripre", "adequac", "dbirmon", "dcntyfipb", "age", "educ")) {
    
    COV <- colnames(smoking) %>% .[grepl(COV, .)]
    
  } 
  
  ### Calculate mean difference effect without covariate
  small_est <- 
    npcausal::ate(y = smoking$dbirwt,
                  a = smoking$is_smoke,
                  x = select(smoking, -dbirwt, -is_smoke) %>% select(-all_of(COV)),
                  nsplits = 5,
                  sl.lib = c("SL.mean", "SL.ranger"))
  
  small_ate_est <- small_est$res$est[3]
  
  ### Calculate measured confounding
  measured_confounding <- abs(full_ate_est - small_ate_est)
  
  cat("Covariate: ", COV, 
      "\nConfounding estimate: ", measured_confounding,
      "\nMaximum: ", max_confounding, "\n")
  
  ### Update covariate corresponding to measured confounding
  if (measured_confounding > max_confounding) {
    max_covariate <- COV
    max_sign <- (full_ate_est - small_ate_est) / abs(full_ate_est - small_ate_est) 
    max_confounding <- measured_confounding
    maxcov_ifvals <- small_est$ifvals$a1 - small_est$ifvals$a0
  }
  
}

confounding_ci_lb <- max_confounding - 
  qnorm(0.975) * sqrt(var(full_ifvals - maxcov_ifvals) / nrow(smoking))

confounding_ci_ub <- max_confounding + 
  qnorm(0.975) * sqrt(var(full_ifvals - maxcov_ifvals) / nrow(smoking))

cat("Covariate corresponding to maximum measured confounding: ", max_covariate,
    "\nMaximum measured confounding: ", max_confounding,
    "\n95% CI: [", confounding_ci_lb, ",", confounding_ci_ub, "]",
    sep = "")

cat("ATE Estimate: ", full_ate_est, 
    "\n95% CI: [", full_ate_lb, full_ate_ub, "]")


##############################################
### 1. Create pointwise confidence interval 
### in sensitivity parameter and find Gamma*
##############################################

########################################
### Traditional sensitivity analysis

dat_old <- data.frame(gamma = seq(from = 0, to = 5 * max_confounding, by = max_confounding / 100))
dat_old$est <- full_ate_est
dat_old$upper_bd_est <- dat_old$est + dat_old$gamma
dat_old$lower_bd_est <- dat_old$est - dat_old$gamma
dat_old$ci_ub <- full_ate_ub + dat_old$gamma
dat_old$ci_lb <- full_ate_lb - dat_old$gamma
dat_old$gamma <- dat_old$gamma / max_confounding
dat_old$analysis <- "old"

Gamma_ast <- min(abs(full_ate_ub), abs(full_ate_lb))

########################################
### Calibrated sensitivity model

dat <- data.frame(gamma = seq(0, 5, 0.01))
dat$est <- full_ate_est
dat$upper_bd_est <- rep(0, nrow(dat))
dat$lower_bd_est <- rep(0, nrow(dat))
dat$ci_lb <- rep(0, nrow(dat))
dat$ci_ub <- rep(0, nrow(dat))

for (row in 1:nrow(dat)) {
  
  dat$upper_bd_est[row] <- full_ate_est + dat$gamma[row] * max_confounding
  dat$lower_bd_est[row] <- full_ate_est - dat$gamma[row] * max_confounding
  
  # When Gamma equals zero, use original analysis
  if (dat$gamma[row] == 0) {
    
    dat$ci_lb[row] <- full_ate_lb
    dat$ci_ub[row] <- full_ate_ub
    
    
    # Otherwise, construct CI based on estimators for upper and lower bounds
  } else {
    
    var_est_ub <- var(full_ifvals + dat$gamma[row] * max_sign * (full_ifvals - maxcov_ifvals))
    var_est_lb <- var(full_ifvals - dat$gamma[row] * max_sign * (full_ifvals - maxcov_ifvals))
    
    dat$ci_lb[row] <- dat$lower_bd_est[row] - qnorm(0.975) * sqrt(var_est_lb / nrow(smoking))
    dat$ci_ub[row] <- dat$upper_bd_est[row] + qnorm(0.975) * sqrt(var_est_ub / nrow(smoking))
    
  }
}

### Calculate rough estimate of Gamma_astM
Gamma_astM <- dat %>% mutate(prod = ci_lb * ci_ub) %>%
  filter(prod == min(abs(prod))) %$% gamma

### Combine with old data and plot bounds
dat$analysis <- "new"
dat %<>% bind_rows(dat_old)

save(dat, smoking, full_ate_est, full_ate_lb, full_ate_ub, full_ifvals, 
     max_confounding, maxcov_ifvals, max_sign, max_covariate, Gamma_ast, Gamma_astM,
     file = "../../intermediate/smoking/smoking_intermediate.RData")


# --------------------------------------------
# 2. Create confidence interval for 
# one-number summary of sensitivity
# --------------------------------------------

onenum_est <- sqrt(full_ate_est^2 / max_confounding^2)
diff_ifvals <- full_ifvals - maxcov_ifvals

onenum_var_est_new <- var(
  (full_ifvals / mean(diff_ifvals)) - (full_ate_est / max_confounding^2) * diff_ifvals
)

onenum_lb_new <- onenum_est - qnorm(0.975) * sqrt(onenum_var_est_new / nrow(smoking))
onenum_ub_new <- onenum_est + qnorm(0.975) * sqrt(onenum_var_est_new / nrow(smoking))

cat("One-number summary of sensitivity: ", onenum_est,
    "\n95% confidence interval: [", onenum_lb_new, ",", onenum_ub_new, "]")

### One-number summary without accounting for measured confounding
onenum_var_est_old <- var(full_ifvals / max_confounding)

onenum_lb_old <- onenum_est - qnorm(0.975) * sqrt(onenum_var_est_old / nrow(smoking))
onenum_ub_old <- onenum_est + qnorm(0.975) * sqrt(onenum_var_est_old / nrow(smoking))

cat("One-number summary of sensitivity (non-ACS): ", onenum_est,
    "\n95% confidence interval: [", onenum_lb_old, ",", onenum_ub_old, "]")

rm(list = ls(all = T))
gc()
