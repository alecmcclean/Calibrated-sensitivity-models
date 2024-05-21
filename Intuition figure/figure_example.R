library(tidyverse)
library(magrittr)
library(latex2exp)
library(ggthemes)

correlation <- seq(-1, 0, 0.01)
rel_std_accuracy <- seq(0, 2, 0.02)
dat <- expand.grid(correlation = correlation, 
                   rel_std_accuracy = rel_std_accuracy)
dat$ratio <- 1 + (dat$rel_std_accuracy)^2 + 2 * dat$correlation * dat$rel_std_accuracy

p <- ggplot(dat = dat) + geom_tile(aes(x = rel_std_accuracy, 
                                  y = correlation, 
                                  color = ratio,
                                  fill = ratio)) +
  scale_fill_gradient2(low = "blue", 
                       mid = "white", 
                       high = "red", 
                       midpoint = 1) +
  scale_color_gradient2(low = "blue", 
                       mid = "white", 
                       high = "red", 
                       midpoint = 1) +
  labs(x = TeX("$\\sqrt{ \\frac{V(\\hat{M}) / M^2}{V \\{ \\hat{u}(\\gamma) \\} / \\gamma^2}}$"),
       y = TeX("$\\rho(\\hat{u}(\\gamma), \\hat{M})$"),
       color = TeX("$\\frac{V \\{ \\hat{U}(\\Gamma) \\}}{V\\{ \\hat{u}(\\gamma) \\}}$"),
       fill =  TeX("$\\frac{V \\{ \\hat{U}(\\Gamma) \\}}{V\\{ \\hat{u}(\\gamma) \\}}$")) +
  theme_clean() +
  theme(legend.position = "top", 
        legend.direction = "horizontal",
        legend.key.height = unit(0.6, "cm"),  # Adjust the height of the legend key
        legend.key.width = unit(2, "cm"),
        axis.text.x = element_text(angle = 90, vjust = 0.2)) +
  geom_text(aes(x = 0, y = -0.8), 
            label = "Sensitivity + post hoc calibration \n under-estimates robustness",
            hjust = 0) +
  geom_text(aes(x = 1.3, y = -0.1),
            label = "Sensitivity + post hoc calibration \n over-estimates robustness",
            hjust = 0)

ggsave(plot = p, filename = "../Figures/intuition.png",
       width = 8, height = 5)

###################################
### Simulations

expit <- function(x) 1 / (1 + exp(-x))
logit <- function(p) log(p / (1-p))

### E(Y^1) = 0.5
### E(E(Y | A = 1, X)) = 0.615525 
### M = 0.213395
### Gamma0 = 1.20887
### Gamma0 / M = 5.7689063422
Gamma0M <- 5.7689063422

N <- 10000
NUM_ITERS <- 100

cis_dat <- data.frame(iter = seq(1:NUM_ITERS))
cis_dat$old_cover <- NA; cis_dat$new_cover <- NA

ests <- data.frame(iter = seq(1:NUM_ITERS))
ests$Gamma0hat <- NA; ests$Mhat <- NA

checker <- data.frame(iter = seq(1:NUM_ITERS))
checker$pt_est <- NA; checker$ci_lb <- NA; checker$ci_ub <- NA

for (ITER in cis_dat$iter) {
  
  dat <- data.frame(X = runif(N, min = -1, max = 1),
                    U = runif(N, min = -1, max = 1))
  dat %<>%
    mutate(
      piXU = (3 + X + U) / 6,
      piX = (3 + X) / 6,
      A = rbinom(N, size = 1, prob = piXU),
      Y1 = 0.5 + X + U,
      Y0 = Y1 - 0.5,
      Y = A * Y1 + (1 - A) * Y0,
      mu1X = 0.5 + X + (1 / (3 * X + 9))
    )
  
  dat$Y1 <- dat$Y1
  dat$Y0 <- dat$Y0
  
  ### "Estimate" relevant nuisance functions
  dat$mu1hat <- dat$mu1X + 0 * N^(-1/2) * rnorm(N)
  
  ### Estimator and variance for M
  if_M <- dat$mu1hat - N * (dat$Y * (dat$A == 1)) / (sum(dat$A == 1)) # Estimator for E(Y | A = 1)
  Mhat <- abs(mean(if_M)); ests$Mhat[ests$iter == ITER] <- Mhat

  ### Estimator and variance for Gamma0:
  Gamma0hat <- abs(mean(dat$mu1hat))
  if_Gamma0 <- dat$mu1hat

  ### Construct confidence intervals for Gamma0 / M, 
  ### with and without incorporating error due to M
  if_old <- if_Gamma0 / Mhat
  if_new <- if_old - (Gamma0hat / Mhat) * ifM
  
  old_cilb <- (Gamma0hat / Mhat) - 1.96 * sqrt(var(if_old) / N) 
  old_ciub <- (Gamma0hat / Mhat) + 1.96 * sqrt(var(if_old) / N)
  
  new_cilb <- (Gamma0hat / Mhat) - 1.96 * sqrt(var(if_new) / N)
  new_ciub <- (Gamma0hat / Mhat) + 1.96 * sqrt(var(if_new) / N)
  
  cis_dat$old_cover[cis_dat$iter == ITER] <- 
    (old_cilb <= Gamma0M) * (Gamma0M <= old_ciub)

  cis_dat$new_cover[cis_dat$iter == ITER] <- 
    (new_cilb <= Gamma0M) * (Gamma0M <= new_ciub)
}

cat("Coverage old: ", mean(cis_dat$old_cover),
    "\nCoverage new: ", mean(cis_dat$new_cover))
