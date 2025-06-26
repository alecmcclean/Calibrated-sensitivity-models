####################################################
### Darfur effect differences analysis
####################################################

sink("../intermediate/darfur/darfur_results.txt")
set.seed(20250624)

###################################
### Load data

load("../intermediate/darfur/darfur_data.RData")

#############################################
### Effect Difference model

### Estimate ATE with full data 
full_est <- 
  npcausal::ate(y = darfur$peacefactor,
                a = darfur$directlyharmed,
                x = darfur_covariates,
                nsplits = 5,
                sl.lib = c("SL.mean", "SL.ranger"))

full_ifvals <- full_est$ifvals$a1 - full_est$ifvals$a0
full_ate_est <- full_est$res$est[3]
full_ate_lb <- full_est$res$ci.ll[3]
full_ate_ub <- full_est$res$ci.ul[3]

# Print results
cat("ATE: ", full_ate_est, 
    "\n95% CI: [", full_ate_lb, ", ", full_ate_ub, "]")


# ---------------------------------------
# Maximum leave-one-out analysis
# ---------------------------------------

### Keep track of covariate corresponding to max confounding
max_confounding <- 0
max_covariate <- "none"

for (COV in c("female", "age", "farmer_dar", "herder_dar", 
              "pastvoted", "hhsize_darfur", "village")) {

  ### Calculate mean difference effect without covariate
  if (COV == "village") {
    COV <- colnames(darfur_covariates) %>% .[grepl("village", .)]
  }
  
  small_est <- 
    npcausal::ate(y = darfur$peacefactor,
                  a = darfur$directlyharmed,
                  x = darfur_covariates %>% dplyr::select(-all_of(COV)),
                  nsplits = 5,
                  sl.lib = c("SL.mean", "SL.ranger"))

  small_ate_est <- small_est$res$est[3]
  small_ifvals   <- small_est$ifvals$a1 - small_est$ifvals$a0
  
  ### Calculate measured confounding
  measured_confounding <- abs(full_ate_est - small_ate_est)
  
  ### SE & 95% CI for this measured confounding
  diff_ifvals    <- full_ifvals - small_ifvals
  se_conf        <- sqrt(var(diff_ifvals) / nrow(darfur))
  ci_lb          <- measured_confounding - qnorm(0.975) * se_conf
  ci_ub          <- measured_confounding + qnorm(0.975) * se_conf
  
  ### Print output
  cat(
    "Covariate:            ", paste(COV, collapse = ","), "\n",
    "Confounding estimate: ", formatC(measured_confounding, digits = 3), "\n",
    "95% CI:               [", formatC(ci_lb, digits = 3),
    ", ", formatC(ci_ub, digits = 3), "]\n\n",
    sep = ""
  )
  
  ### Update covariate corresponding to measured confounding
  if (measured_confounding > max_confounding) {
    max_covariate <- COV
    max_sign <- (full_ate_est - small_ate_est) / abs(full_ate_est - small_ate_est) 
    max_confounding <- measured_confounding
    maxcov_ifvals <- small_est$ifvals$a1 - small_est$ifvals$a0
  }

}

confounding_ci_lb <- max_confounding - 
  qnorm(0.975) * sqrt(var(full_ifvals - maxcov_ifvals) / nrow(darfur))

confounding_ci_ub <- max_confounding + 
  qnorm(0.975) * sqrt(var(full_ifvals - maxcov_ifvals) / nrow(darfur))

cat("Covariate corresponding to maximum measured confounding: ", max_covariate,
    "\nMaximum measured confounding: ", max_confounding,
    "\n95% CI: [", confounding_ci_lb, ",", confounding_ci_ub, "]",
    sep = "")
    
cat("ATE Estimate: ", full_ate_est, 
    "\n95% CI: [", full_ate_lb, full_ate_ub, "]")


# --------------------------------------------
# 1. Create pointwise confidence interval 
# in sensitivity parameter and find Gamma*
# --------------------------------------------

# --------------------------------------------
# Traditional sensitivity analysis

dat_old <- data.frame(gamma = seq(from = 0, to = 5 * max_confounding, by = max_confounding / 100))
dat_old$est <- full_ate_est
dat_old$upper_bd_est <- dat_old$est + dat_old$gamma
dat_old$lower_bd_est <- dat_old$est - dat_old$gamma
dat_old$ci_ub <- full_ate_ub + dat_old$gamma
dat_old$ci_lb <- full_ate_lb - dat_old$gamma
dat_old$gamma <- dat_old$gamma / max_confounding
dat_old$analysis <- "old"

Gamma_ast <- min(abs(full_ate_ub), abs(full_ate_lb))

# --------------------------------------------
# ACS

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
  # See Theorem 1
  } else {
    
    var_est_ub <- var(full_ifvals + dat$gamma[row] * max_sign * (full_ifvals - maxcov_ifvals))
    var_est_lb <- var(full_ifvals - dat$gamma[row] * max_sign * (full_ifvals - maxcov_ifvals))
    
    dat$ci_ub[row] <- dat$upper_bd_est[row] + qnorm(0.975) * sqrt(var_est_ub / nrow(darfur))
    dat$ci_lb[row] <- dat$lower_bd_est[row] - qnorm(0.975) * sqrt(var_est_lb / nrow(darfur))
    
  }
}


### Calculate rough estimate of Gamma_astM
Gamma_astM <- dat %>% mutate(prod = ci_lb * ci_ub) %>%
  filter(prod == min(abs(prod))) %$% gamma

### Combine with old data and plot bounds
dat$analysis <- "new"
dat %<>% bind_rows(dat_old)

save(dat, darfur, full_ate_est, full_ate_lb, full_ate_ub, full_ifvals, 
     max_confounding, maxcov_ifvals, max_sign, max_covariate, Gamma_ast, Gamma_astM,
     file = "../intermediate/darfur/fx_intermediate.RData")


# --------------------------------------------
# Print results for point estimates
# --------------------------------------------

# \Gamma_\ast / M from standard approach

cat("Gamma_ast / M from usual approach: ", Gamma_ast / max_confounding,
    "\nGamma_ast^M from new approach: ", Gamma_astM)

# --------------------------------------------
# 2. Create confidence interval for 
# one-number summary of sensitivity
# --------------------------------------------

onenum_est <- sqrt(full_ate_est^2 / max_confounding^2)
diff_ifvals <- full_ifvals - maxcov_ifvals

onenum_var_est_new <- var(
  (full_ifvals / mean(diff_ifvals)) - (full_ate_est / max_confounding^2) * diff_ifvals
)

onenum_lb_new <- onenum_est - qnorm(0.975) * sqrt(onenum_var_est_new / nrow(darfur))
onenum_ub_new <- onenum_est + qnorm(0.975) * sqrt(onenum_var_est_new / nrow(darfur))

cat("One-number summary of sensitivity: ", onenum_est,
    "\n95% confidence interval: [", onenum_lb_new, ",", onenum_ub_new, "]")

### One-number summary without accounting for measured confounding
onenum_var_est_old <- var(full_ifvals / max_confounding)

onenum_lb_old <- onenum_est - qnorm(0.975) * sqrt(onenum_var_est_old / nrow(darfur))
onenum_ub_old <- onenum_est + qnorm(0.975) * sqrt(onenum_var_est_old / nrow(darfur))

cat("One-number summary of sensitivity (non-ACS): ", onenum_est,
    "\n95% confidence interval: [", onenum_lb_old, ",", onenum_ub_old, "]")

save(onenum_est, onenum_var_est_new, onenum_lb_new, onenum_ub_new,
     onenum_var_est_old, onenum_lb_old, onenum_ub_old, 
     Gamma_ast, max_confounding, Gamma_astM, file = "../data/darfur_results.RData")

sink()