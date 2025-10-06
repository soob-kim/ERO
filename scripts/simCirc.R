# ==========================================================
# File: simCirc.R
# Purpose: Simulation in circular data (Section 3.1 Fig 3)
# Author: Soobin Kim
# ==========================================================
# - Sample size: n = 100
# - Distribution:
#     theta = seq(0, 2pi), rho \sim N(2, 0.1^2)
# ---------------------------------------------------------
source("R/euclidean.R")

library(dplyr)
library(ggplot2)
library(reshape2)
library(gridExtra)
library(plotly)

# Generate data
set.seed(150)
n = 100
theta = seq(0, 2*pi, length = n)
rho = rnorm(n, mean = 1, sd = 0.05)
simcirc = cbind(2 * rho * cos(theta), 2 * rho * sin(theta))

# Calculate ball statistics
simcirc_BMT <- EucBMTraj(simcirc)
simcirc_DT <- EucDistTraj(simcirc_BMT)

# FPCA of distance trajectory
simcirc_DFPCA <- EucDFPCA(simcirc_DT, ind = factor(1:n),
                          legend_loc = c(0.23, 0.23),
                          yliml = -0.01, ylimu = 0.5)

# Plot
#### Figure 3
fig3_1 <- data.frame(simcirc) %>%
  ggplot(data = ., aes(x = X1, y = X2)) +
  geom_point() +
  theme_light()

fig3_2 <- plot_ly(type = 'scatter3d') 
fig3_2 <- add_trace(
  fig3_2, mode = "lines",
  x = simcirc_BMT$traj[, "r"], 
  y = simcirc_BMT$traj[, "X1"], 
  z = simcirc_BMT$traj[, "X2"],
  color = as.factor(simcirc_BMT$traj[, "x0"]),
  line = list(color = "gray", opacity = 1, width = 0.7), showlegend = F)
fig3_2 <- add_trace(
  fig3_2, mode = "markers", x = rep(0, n), 
  y = simcirc[,1],
  z = simcirc[,2],
  color = as.factor(simcirc_BMT$traj[which(simcirc_BMT$traj[, "r"] == 0), "x0"]),
  marker = list(size = 3, color = "midnightblue"), showlegend = F) %>%
  layout(scene = list(xaxis = list(title = 'r'),
                      yaxis = list(title = 'X1'),
                      zaxis = list(title = 'X2'),
                      aspectratio = list(x=1, y=1, z=1)))



fig3_3 <- as_tibble(simcirc_DT$traj, 
                    .name_repair = ~ as.character(simcirc_DT$grid)) %>%
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


fig3_6 <- simcirc_DFPCA$fpc %>%
  ggplot(data = ., aes(x = X1, y = X2)) +
  geom_point() +
  xlab("PC1") + ylab("PC2") +
  guides(color = "none") +
  theme_light()


fig3_1
fig3_2
fig3_3
simcirc_DFPCA$fig_mov1
simcirc_DFPCA$fig_mov2
fig3_6
