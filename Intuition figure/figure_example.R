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
  labs(x = TeX("$\\frac{\\sqrt{V(\\hat{M})} / M}{\\sqrt{V \\{ \\hat{u}(\\gamma) \\}} / \\gamma}$"),
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

ggsave(plot = p, filename = "../figures/intuition.png",
       width = 8, height = 5)
