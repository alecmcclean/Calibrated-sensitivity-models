####################################################
### Figures for smoking data analysis
####################################################

load("../intermediate/smoking/smoking_intermediate.RData")

options(ggplot2.discrete.colour= c("red", "blue"))

dat %<>% mutate(
  analysis = ifelse(analysis == "new",
                    "Calibrated sensitivity",
                    "Standard sensitivity")
)

### Plot bounds
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

ggsave("../figures/smoking_birthweight_fx.png", plot = p,
       width = 8, height = 6)

options(ggplot2.discrete.colour= c("blue", "red"))

p <- ggplot(data = dat %>% filter(analysis == "Standard sensitivity"), 
            aes(x = gamma * max_confounding)) + 
  geom_line(aes(y = est), linewidth = 1, linetype = "dashed") +
  geom_line(aes(y = upper_bd_est, color = analysis), color = "purple", linewidth = 1) +
  geom_line(aes(y = lower_bd_est, color = analysis), color = "purple", linewidth = 1) +
  geom_line(aes(y = ci_lb, color = analysis), linewidth = 1, linetype = "dashed") +
  geom_line(aes(y = ci_ub, color = analysis), linewidth = 1, linetype = "dashed") +
  geom_ribbon(aes(ymin = ci_lb, ymax = ci_ub, color = analysis), alpha = 0.2) +
  geom_hline(yintercept = 0, linewidth = 1, linetype = "dotted") +
  geom_vline(xintercept = max_confounding, linewidth = 1, linetype = "dashed", color = "blue") +
      theme_clean(base_size = 15) +
  scale_y_continuous(limits = c(-1000, 250), breaks = seq(-1000, 250, by = 250)) +
  geom_text(aes(x = max_confounding - 13, y = -500), 
            label = TeX("$\\hat{M} = 46.6$"), 
            color = "blue") +
  theme(legend.position = "top", 
        legend.direction = "horizontal",
        legend.key.height = unit(0.6, "cm"),  # Adjust the height of the legend key
        legend.key.width = unit(1.5, "cm"),
        legend.text = element_text(size = 12),
        axis.text.y = element_text(angle = 60)) +
  labs(x = TeX("$\\gamma$ in sensitivity model"),
       y = TeX("ATE, bounds, and pointwise 95\\% CI"),
       color = "Analysis")

ggsave("../figures/smoking_intro_standard_sensitivity.png", plot = p,
       width = 8, height = 6)

options(ggplot2.discrete.colour= c("red"))

p <- ggplot(data = dat %>% filter(analysis == "Calibrated sensitivity"), 
            aes(x = gamma)) + 
  geom_line(aes(y = est), linewidth = 1, linetype = "dashed") +
  geom_line(aes(y = upper_bd_est, color = analysis), color = "purple", linewidth = 1) +
  geom_line(aes(y = lower_bd_est, color = analysis), color = "purple", linewidth = 1) +
  geom_line(aes(y = ci_lb, color = analysis), linewidth = 1, linetype = "dashed") +
  geom_line(aes(y = ci_ub, color = analysis), linewidth = 1, linetype = "dashed") +
  geom_ribbon(aes(ymin = ci_lb, ymax = ci_ub, color = analysis), alpha = 0.2) +
  geom_hline(yintercept = 0, linewidth = 1, linetype = "dotted") +
  scale_y_continuous(limits = c(-1000, 250), breaks = seq(-1000, 250, by = 250)) +
  theme_clean(base_size = 15) +
  theme(legend.position = "top", 
        legend.direction = "horizontal",
        legend.key.height = unit(0.6, "cm"),  # Adjust the height of the legend key
        legend.key.width = unit(1.5, "cm"),
        legend.text = element_text(size = 12),
        axis.text.y = element_text(angle = 60)) +
  labs(x = TeX("$\\Gamma$ in calibrated sensitivity model"),
       y = TeX("ATE, bounds, and pointwise 95\\% CI"),
       color = "Analysis")

ggsave("../figures/smoking_intro_calibrated_sensitivity.png", plot = p,
       width = 8, height = 6)

rm(list = ls(all = T))
gc()
