####################################################
### Figures for Darfur data analysis
####################################################

load("../intermediate/darfur/fx_intermediate.RData")

ATE_EST <- dat$est[1]
dat %<>% select(-est) %>% 
  gather(var, value, upper_bd_est, lower_bd_est, ci_lb, ci_ub) %>% 
  filter(!(analysis == "new" & var == "lower_bd_est"),
         !(analysis == "new" & var == "upper_bd_est")) %>% 
  mutate(analysis = ifelse(grepl("_est", var), "bound", analysis),
         var = ifelse(var %in% c("ci_lb", "lower_bd_est"), "lower", "upper")) 

dat %<>% mutate(
  analysis = case_when(
    analysis == "new" ~ "Calibrated sensitivity CI",
    analysis == "old" ~ "Standard sensitivity CI",
    T ~ "Estimated bounds"
  ),
  analysis = factor(analysis,
                    levels = c("Estimated bounds", 
                               "Standard sensitivity CI",
                               "Calibrated sensitivity CI"))
)

dat %<>% spread(var, value)

################################
### Effect differences

p <- ggplot(data = dat, aes(x = gamma)) + 
  geom_line(aes(y = ATE_EST), linewidth = 1, linetype = "dashed") +
  geom_line(aes(y = upper, color = analysis, linetype = analysis), linewidth = 1) +
  geom_line(aes(y = lower, color = analysis, linetype = analysis), linewidth = 1) +
  geom_hline(yintercept = 0, linewidth = 1, linetype = "dotted") +
  theme_clean(base_size = 15) +
  scale_color_manual(values = c("purple", "blue", "red")) +
  scale_linetype_manual(values = c("dotdash", "longdash", "solid")) +
  scale_x_continuous(breaks = seq(0, 5, 1)) +
  scale_y_continuous(limits = c(-0.12, 0.22), n.breaks = 10) +
  theme(legend.position = "top", 
        legend.direction = "horizontal",
        legend.key.height = unit(0.6, "cm"),  # Adjust the height of the legend key
        legend.key.width = unit(1.5, "cm"),
        legend.text = element_text(size = 10),
        axis.text.y = element_text(angle = 60)) +
  labs(x = TeX("$\\Gamma$ in calibrated sensitivity model, $\\frac{\\gamma}{\\hat{M}}$ in sensitivity model"),
       y = TeX("ATE, bounds, and pointwise 95\\% CI"),
       color = "", 
       linetype = "")

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

plot_dat %<>%
  select(analysis, gamma = Gamma, upper_bd_est = upper_est, 
         lower_bd_est = lower_est, ci_ub = upper_ci_ub, ci_lb = lower_ci_lb)
  
plot_dat %<>% 
  gather(var, value, upper_bd_est, lower_bd_est, ci_lb, ci_ub) %>% 
  filter(!(analysis == "new" & var == "lower_bd_est"),
         !(analysis == "new" & var == "upper_bd_est")) %>% 
  mutate(analysis = ifelse(grepl("_est", var), "bound", analysis),
         var = ifelse(var %in% c("ci_lb", "lower_bd_est"), "lower", "upper")) 

plot_dat %<>% mutate(
  analysis = case_when(
    analysis == "new" ~ "Calibrated sensitivity CI",
    analysis == "old" ~ "Standard sensitivity CI",
    T ~ "Estimated bounds"
  ),
  analysis = factor(analysis,
                    levels = c("Estimated bounds", 
                               "Standard sensitivity CI",
                               "Calibrated sensitivity CI"))
)

plot_dat %<>% spread(var, value)

########################################
### Plots!

ATE_EST <- full_ate_est

p <- ggplot(data = plot_dat, aes(x = log(gamma) / estimates$max_confounding[1])) + 
  geom_line(aes(y = ATE_EST), linewidth = 1, linetype = "dashed") +
  geom_line(aes(y = upper, color = analysis, linetype = analysis), linewidth = 1) +
  geom_line(aes(y = lower, color = analysis, linetype = analysis), linewidth = 1) +
  geom_hline(yintercept = 0, linewidth = 1, linetype = "dotted") +
  theme_clean(base_size = 15) +
  scale_color_manual(values = c("purple", "blue", "red")) +
  scale_linetype_manual(values = c("dotdash", "longdash", "solid")) +
  scale_y_continuous(limits = c(-0.12, 0.22), n.breaks = 10) +
  theme(legend.position = "top", 
        legend.direction = "horizontal",
        legend.key.height = unit(0.6, "cm"),  # Adjust the height of the legend key
        legend.key.width = unit(1.5, "cm"),
        legend.text = element_text(size = 10),
        axis.text.y = element_text(angle = 60)) +
  labs(x = TeX("$\\Gamma$ in calibrated sensitivity model, $\\frac{\\gamma}{\\hat{M}}$ in sensitivity model"),
       y = TeX("ATE, bounds, and pointwise 95\\% CI"),
       color = "", 
       linetype = "")

ggsave(plot = p, filename = "../figures/violence_peace_odds.png",
       width = 8, height = 6)

rm(list = ls(all = T))
gc()
