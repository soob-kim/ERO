# ==========================================================
# File: simR2.R
# Purpose: Simulation in R^2 (Section 3.1 Fig 1, 2)
# Author: Soobin Kim
# ==========================================================
# - Sample size: n = 100
# - Distribution:
#     (X1, X2) \sim N((0,0), diag(c(1, 1/2)))
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
n <- 100
simr2 <- MASS::mvrnorm(n, mu = c(0,0), Sigma = diag(c(1, 1/2)))

# Calculate ball statistics
simr2_BMT <- EucBMTraj(simr2)
simr2_DT <- EucDistTraj(simr2_BMT)

# FPCA of distance trajectory
simr2_DFPCA <- EucDFPCA(simr2_DT, ind = factor(1:n))

# Plot
#### Figure 1
fig1 <- plot_ly(type = 'scatter3d') 
fig1 <- add_trace(
  fig1, mode = "lines",
  x = simr2_BMT$traj[, "r"], 
  y = simr2_BMT$traj[, "X1"], 
  z = simr2_BMT$traj[, "X2"],
  color = as.factor(simr2_BMT$traj[, "x0"]),
  line = list(color = "gray", opacity = 1, width = 0.7), showlegend = F)
fig1 <- add_trace(
  fig1, mode = "markers", x = rep(0, n), 
  y = simr2[,1],
  z = simr2[,2],
  color = as.factor(simr2_BMT$traj[which(simr2_BMT$traj[, "r"] == 0), "x0"]),
  marker = list(size = 3, color = "midnightblue"), showlegend = F) %>%
  layout(scene = list(xaxis = list(title = 'r'),
                      yaxis = list(title = 'X1'),
                      zaxis = list(title = 'X2')))

#### Figure 2
fig2_1 <- as_tibble(simr2_DT$traj, 
          .name_repair = ~ as.character(simr2_DT$grid)) %>%
  mutate(ind = 1:n) %>%
  melt(id.vars = "ind", value.name = "d", variable.name = "r") %>%
  mutate(r = as.numeric(levels(r))[r], d = d,
         ind = as.numeric(ind)) %>%
  ggplot(aes(x = r, y = d, group = ind)) +
  ylim(0, 0.6) +
  geom_line(alpha = 0.5) +
  scale_x_continuous(breaks = seq(0, 1, by = 0.2),
                     limits = c(0, 1.1)) +
  theme_light()

fig2_3 <- simr2_DFPCA$fpc %>%
  ggplot(data = ., aes(x = X1, y = X2)) +
  geom_point() +
  xlab("PC1") + ylab("PC2") +
  guides(color = "none") +
  theme_light()

grid.arrange(fig2_1, simr2_DFPCA$fig_mov1, fig2_3, simr2_DFPCA$fig_mov2, nrow = 2)

