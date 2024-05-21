####################################################
### Figures for Darfur data analysis
####################################################

load("../intermediate/darfur/fx_intermediate.RData")

dat %<>% mutate(
  analysis = ifelse(analysis == "new",
                    "Calibrated sensitivity",
                    "Standard sensitivity")
)

################################
### Effect differences

options(ggplot2.discrete.colour= c("red", "blue"))
p <- ggplot(data = dat, aes(x = gamma)) + 
  geom_line(aes(y = est), linewidth = 1, linetype = "dashed") +
  geom_line(aes(y = upper_bd_est, color = analysis), color = "purple", linewidth = 1) +
  geom_line(aes(y = lower_bd_est, color = analysis), color = "purple", linewidth = 1) +
  geom_line(aes(y = ci_lb, color = analysis), linewidth = 1, linetype = "dashed") +
  geom_line(aes(y = ci_ub, color = analysis), linewidth = 1, linetype = "dashed") +
  geom_ribbon(aes(ymin = ci_lb, ymax = ci_ub, color = analysis), alpha = 0.2) +
  geom_hline(yintercept = 0, linewidth = 1, linetype = "dotted") +
  theme_clean(base_size = 15) +
  scale_x_continuous(breaks = seq(0, 5, 1)) +
  theme(legend.position = "top", 
        legend.direction = "horizontal",
        legend.key.height = unit(0.6, "cm"),  # Adjust the height of the legend key
        legend.key.width = unit(1.5, "cm"),
        legend.text = element_text(size = 12),
        axis.text.y = element_text(angle = 60)) +
  labs(x = TeX("$\\Gamma$ in calibrated sensitivity model, $\\frac{\\gamma}{\\hat{M}}$ in sensitivity model"), 
       y = TeX("ATE, bounds, and pointwise 95\\% CI"),
       color = "Analysis")

ggsave(plot = p, filename = "../figures/violence_peace_fx.png",
       width = 8, height = 6)

rm(list = ls(all = T))
gc()

##############################
### Bounded odds ratio

estimates <- readRDS("../intermediate/darfur/bdd_odds_no_M.RDS")
bootstraps <- readRDS("../intermediate/darfur/bdd_odds_bootstraps.RDS")

###########################
### Measured confounding

M_sd <- bootstraps %>% 
  filter(Gamma == 1) %>% 
  select(max_confounding) %$% 
  sd(max_confounding)

cat("Measured confounding: ", unique(estimates$max_confounding),
    " in [", unique(estimates$max_confounding) - qnorm(0.975) * M_sd,
    ", ", unique(estimates$max_confounding) + qnorm(0.975) * M_sd, "]")

### Compute variance estimate from bootstrap
se_estimates <- bootstraps %>% 
  gather(var, value, upper_est, lower_est) %>% 
  group_by(Gamma, var) %>% 
  summarize(se_est = sqrt(var(value))) %>%
  spread(var, se_est)

# Construct confidence intervals including uncertainty in M
for (gamma in unique(estimates$Gamma)) {
  
  if (gamma == 1) {
    estimates$upper_ci_M[estimates$Gamma == gamma] <-
      estimates$upper_ci_ub[estimates$Gamma == gamma] 
    
    estimates$lower_ci_M[estimates$Gamma == gamma] <-
      estimates$lower_ci_lb[estimates$Gamma == gamma] 

  } else {
  
    estimates$upper_ci_M[estimates$Gamma == gamma] <-
      estimates$upper_est[estimates$Gamma == gamma] + 
      qnorm(0.975) * se_estimates$upper_est[se_estimates$Gamma == gamma]
  
    estimates$lower_ci_M[estimates$Gamma == gamma] <-
      estimates$lower_est[estimates$Gamma == gamma] - 
      qnorm(0.975) * se_estimates$lower_est[se_estimates$Gamma == gamma] 
  }
}

plot_dat <- estimates %>% select(-upper_ci_ub, -lower_ci_lb) %>%
  mutate(analysis = "new") %>%
  rename(lower_ci_lb = lower_ci_M, upper_ci_ub = upper_ci_M) %>%
  bind_rows(
    estimates %>% mutate(analysis = "old") %>% select(-upper_ci_M, -lower_ci_M)
  ) %>%
  select(analysis, everything()) %>%
  arrange(analysis, Gamma)

load("../intermediate/darfur/fx_intermediate.RData")
plot_dat %<>% mutate(
  upper_est = ifelse(Gamma == 1, full_ate_est, upper_est),
  lower_est = ifelse(Gamma == 1, full_ate_est, lower_est),
  upper_ci_ub = ifelse(Gamma == 1, full_ate_ub, upper_ci_ub),
  lower_ci_ub = ifelse(Gamma == 1, full_ate_lb, lower_ci_lb)
)
  

########################################
### Plots!

ATE_EST <- plot_dat$upper_est[1]
options(ggplot2.discrete.colour= c("red", "blue"))

# Guesstimate Gamma_+ for now; later, interpolate with model
Gamma_astM <- 0.025
Gamma_ast <- 0.034

p <- plot_dat %>%
  mutate(analysis = ifelse(analysis == "new", "Calibrated sensitivity", "Standard sensitivity")) %>%
  ggplot(aes(x = log(Gamma) / plot_dat$max_confounding[1])) +
  geom_hline(yintercept = ATE_EST, linewidth = 1, linetype = "dashed") +
  geom_line(aes(y = upper_est, color = analysis), color = "purple", linewidth = 1) +
  geom_line(aes(y = lower_est, color = analysis), color = "purple", linewidth = 1) +
  geom_line(aes(y = lower_ci_lb, color = analysis), linewidth = 1, linetype = "dashed") +
  geom_line(aes(y = upper_ci_ub, color = analysis), linewidth = 1, linetype = "dashed") +
  geom_ribbon(aes(ymin = lower_ci_lb, ymax = upper_ci_ub, color = analysis), alpha = 0.2) +
  geom_hline(yintercept = 0, linewidth = 1, linetype = "dotted") +

  theme_clean(base_size = 15) +
  theme(legend.position = "top", 
        legend.direction = "horizontal",
        legend.key.height = unit(0.6, "cm"),  # Adjust the height of the legend key
        legend.key.width = unit(1.5, "cm"),
        legend.text = element_text(size = 12),
        axis.text.y = element_text(angle = 60)) +
  scale_y_continuous(limits = c(-0.15, 0.25), n.breaks = 10) +
  labs(x = TeX("$\\Gamma$ in calibrated sensitivity model, $\\frac{\\gamma}{\\widehat{M}}$ in sensitivity model"),
       y = TeX("ATE, bounds, and pointwise 95\\% CI"),
       color = "Analysis")

ggsave(plot = p, filename = "../figures/violence_peace_odds.png",
       width = 8, height = 6)

rm(list = ls(all = T))
gc()
