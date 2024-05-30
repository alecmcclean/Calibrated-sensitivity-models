####################################################
### Figures for smoking data analysis
####################################################

load("../intermediate/smoking/smoking_intermediate.RData")

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

### Plot bounds
p <- ggplot(data = dat, aes(x = gamma)) + 
  geom_line(aes(y = ATE_EST), linewidth = 1, linetype = "dashed") +
  geom_line(aes(y = upper, color = analysis, linetype = analysis), linewidth = 1) +
  geom_line(aes(y = lower, color = analysis, linetype = analysis), linewidth = 1) +
  geom_hline(yintercept = 0, linewidth = 1, linetype = "dotted") +
  theme_clean(base_size = 15) +
  scale_color_manual(values = c("purple", "blue", "red")) +
  scale_linetype_manual(values = c("dotdash", "longdash", "solid")) +
  scale_x_continuous(breaks = seq(0, 5, 1)) +
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

ggsave("../figures/smoking_birthweight_fx.png", plot = p,
       width = 8, height = 6)

rm(list = ls(all = T))
gc()
