set.seed(20231205)

###########################################
### Install npcausal for ATE estimation

library(npcausal)

###########################################
### Install package to load birth weight data

library(foreign)

############################################################
### Install packages for data manipulation and plots

library(tidyverse)
library(ggthemes) # I like clean_theme
library(latex2exp)


########################
### Load data

smoking <- read.dta("C_2010_JOE-dataRS5K.dta")

### Dichotomize treatment
smoking$is_smoke <- ifelse(smoking$T > 0, 1, 0)
smoking <- select(smoking, -T) 

### Remove constant
smoking <- select(smoking, -const)

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
          "age", "educ", "nprevisit", "disllb",  "dlivord")

for (COV in COVS) {
  
  # Data cleaning --- remove all reasonable groups of variables as if they
  # were one variable (e.g., race, age, education of parents)
  if (COV == "race") {# Remove all race 
    
    COV <- c("mwhite", "mblack", "mhispan", "fwhite", "fblack", "fhispan")
    
  } else if (COV %in% c("age", "educ")) {
    
    COV <- colnames(smoking) %>% .[grepl(COV, .)]
    
  } 
  
  # Remove all binary variables that are levels of the same factor variables together
  if (COV %in% c("tripre", "adequac", "dbirmon", "dcntyfipb", "age", "educ")) {
    
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
### in sensitivity parameter
##############################################

dat <- data.frame(gamma = seq(0, 5, 0.5))
dat$est <- full_ate_est
dat$upper_bd_est <- rep(0, nrow(dat))
dat$lower_bd_est <- rep(0, nrow(dat))
dat$ci_lb <- rep(0, nrow(dat))
dat$ci_ub <- rep(0, nrow(dat))

for (row in 1:nrow(dat)) {
  
  # When Gamma equals zero, use original analysis
  if (dat$gamma[row] == 0) {
    
    dat$ci_lb[row] <- full_ate_lb
    dat$ci_ub[row] <- full_ate_ub
    
    
    # Otherwise, construct CI based on estimators for upper and lower bounds
    # See Theorem 1 in McClean et al. 2023
  } else {
    
    dat$upper_bd_est[row] <- full_ate_est + dat$gamma[row] * max_confounding
    dat$lower_bd_est[row] <- full_ate_est - dat$gamma[row] * max_confounding
    var_est_ub <- var(full_ifvals + dat$gamma[row] * max_sign * (full_ifvals - maxcov_ifvals))
    var_est_lb <- var(full_ifvals - dat$gamma[row] * max_sign * (full_ifvals - maxcov_ifvals))
    
    dat$ci_lb[row] <- dat$lower_bd_est[row] - qnorm(0.975) * sqrt(var_est_lb / nrow(smoking))
    dat$ci_ub[row] <- dat$upper_bd_est[row] + qnorm(0.975) * sqrt(var_est_ub / nrow(smoking))
    
  }
}

p <- ggplot(data = dat, aes(x = gamma, y = est)) + 
  geom_line(size = 1) +
  geom_line(aes(x = gamma, y = ci_lb), size = 1.5, color = "blue") +
  geom_line(aes(x = gamma, y = ci_ub), size = 1.5, color = "blue") +
  geom_ribbon(aes(ymin = ci_lb, ymax = ci_ub), alpha = 0.2) +
  geom_hline(yintercept = 0, size = 1.5, color = "red", linetype = "dashed") +
  theme_clean(base_size = 15) +
  scale_x_continuous(breaks = seq(0, 10, 1)) +
  labs(x = TeX("$\\Gamma$ - sensitivity parameter"), 
       y = TeX("ATE and pointwise 95\\% CI"))

ggsave("data/smoking_birthweight.png", plot = p,
       width = 6, height = 4)

#################################################
### 2. Create confidence interval for 
### one-number summary of sensitivity
################################################

onenum_est <- sqrt(full_ate_est^2 / max_confounding^2)
diff_ifvals <- full_ifvals - maxcov_ifvals

onenum_var_est <- var(
  (full_ifvals / mean(diff_ifvals)) - (full_ate_est / max_confounding^2) * diff_ifvals
)

onenum_lb <- onenum_est - qnorm(0.975) * sqrt(onenum_var_est / nrow(smoking))
onenum_ub <- onenum_est + qnorm(0.975) * sqrt(onenum_var_est / nrow(smoking))

cat("One-number summary of sensitivity: ", onenum_est,
    "\n95% confidence interval: [", onenum_lb, ",", onenum_ub, "]")


