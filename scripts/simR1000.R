# ==========================================================
# File: simR1000.R
# Purpose: Simulation in R^{1000} (Section 3.1 Fig 5)
# Author: Soobin Kim
# ==========================================================
# - Sample size: n = 50
# - Distribution:
#     (X1, .., X1000) \sim N_{1000}((0,...,0), diag(c(1, 1/2,..., 1/1000)))
# ---------------------------------------------------------
source("R/euclidean.R")

library(dplyr)
library(ggplot2)
library(reshape2)
library(fdapace) # for FPCA
library(gridExtra)
library(plotly)

# Generate data
set.seed(100)
n = 50
simhd <- MASS::mvrnorm(n, mu=rep(0, 1000), Sigma=diag(c(1:1000)^(-1)))

# Calculate ball statistics
simhd_BMT <- EucBMTraj(simhd)
simhd_DT <- EucDistTraj(simhd_BMT)
simhd_DT_trunc <- list(traj = simhd_DT$traj[, 21:56],
                       grid = seq(0.4, 1.1, length = 36))
# FPCA of distance trajectory
simhd_DFPCA <- EucDFPCA(simhd_DT_trunc, ind = factor(1:n),
                        opts = NULL,
                        yliml = -0.01, ylimu = 0.55)

# Plot
#### Figure 5
fig5_1 <- as_tibble(simhd_DT_trunc$traj, 
                    .name_repair = ~ as.character(simhd_DT_trunc$grid)) %>%
  mutate(ind = 1:n) %>%
  melt(id.vars = "ind", value.name = "d", variable.name = "r") %>%
  mutate(r = as.numeric(levels(r))[r], d = d,
         ind = as.numeric(ind)) %>%
  ggplot(aes(x = r, y = d, group = ind)) +
  geom_line(alpha = 0.5) +
  scale_x_continuous(breaks = seq(0.4, 1, by = 0.2),
                     limits = c(0.4, 1.1)) +
  theme_light()

fig5_3 <- simhd_DFPCA$fpc %>%
  ggplot(data = ., aes(x = X1, y = X2)) +
  geom_point() +
  xlab("PC1") + ylab("PC2") +
  guides(color = "none") +
  theme_light()

grid.arrange(fig5_1, simhd_DFPCA$fig_mov1, fig5_3, simhd_DFPCA$fig_mov2, nrow = 2)

