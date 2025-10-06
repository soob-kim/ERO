# ==========================================================
# File: simCirc2.R
# Purpose: Simulation in two circular data (Section 3.1 Fig 4)
# Author: Soobin Kim
# ==========================================================
# - Sample size: n = 100
# - Distribution:
#     Group 1: theta = seq(0, 2pi), rho \sim Unif(0,1)
#     Group 2: theta = seq(0, 2pi), rho \sim N(2,0.2^2)
# ---------------------------------------------------------
source("R/euclidean.R")

library(dplyr)
library(ggplot2)
library(reshape2)
library(gridExtra)
library(plotly)

# Generate data
set.seed(100)
n = 100
rad = runif(n/2, 0, 1)
theta = seq(0, 2*pi, length = n/2)
rho = rnorm(n/2, mean = 1, sd = 0.1)
simcirc2 = rbind(
  cbind(rad * cos(theta), rad * sin(theta)),
  cbind(2 * rho * cos(theta), 2 * rho * sin(theta))
)

# Calculate ball statistics
simcirc2_BMT <- EucBMTraj(simcirc2)
simcirc2_DT <- EucDistTraj(simcirc2_BMT)

# FPCA of distance trajectory
simcirc2_DFPCA <- EucDFPCA(simcirc2_DT, ind = factor(1:n),
                          yliml = -0.01, ylimu = 0.45)

# Plot
#### Figure 4
fig4_1 <- data.frame(simcirc2, gr = factor(rep(c(1, 2), each = n/2))) %>%
  ggplot(data = ., aes(x = X1, y = X2, color = gr, shape = gr)) +
  geom_point() +
  scale_color_manual(values=c("black", "tomato")) +
  guides(color = "none", shape = "none") +
  theme_light()

fig4_2 <- plot_ly(type = 'scatter3d') 
fig4_2 <- add_trace(
  fig4_2, mode = "lines",
  x = simcirc2_BMT$traj[, "r"], 
  y = simcirc2_BMT$traj[, "X1"], 
  z = simcirc2_BMT$traj[, "X2"],
  color = as.factor(simcirc2_BMT$traj[, "x0"]),
  line = list(color = "gray", opacity = 1, width = 0.7), showlegend = F)
fig4_2 <- add_trace(
  fig4_2, mode = "markers", x = rep(0, n), 
  y = simcirc2[,1],
  z = simcirc2[,2],
  color = as.factor(simcirc2_BMT$traj[which(simcirc2_BMT$traj[, "r"] == 0), "x0"]),
  marker = list(size = 3, color = c(rep("midnightblue", 50), rep("tomato", 50))), 
  showlegend = F) %>%
  layout(scene = list(xaxis = list(title = 'r'),
                      yaxis = list(title = 'X1'),
                      zaxis = list(title = 'X2'),
                      aspectratio = list(x=1, y=1, z=1)))



fig4_3 <- as_tibble(simcirc2_DT$traj, 
                    .name_repair = ~ as.character(simcirc2_DT$grid)) %>%
  mutate(ind = 1:n) %>%
  melt(id.vars = "ind", value.name = "d", variable.name = "r") %>%
  mutate(r = as.numeric(levels(r))[r], d = d,
         ind = as.numeric(ind),
         gr = factor(if_else(ind <= 50, 1, 2))) %>%
  ggplot(aes(x = r, y = d, group = ind, colour = gr, linetype = gr)) +
  ylim(0, 0.6) +
  geom_line(alpha = 0.5) +
  scale_color_manual(values=c("black", "tomato")) +
  scale_x_continuous(breaks = seq(0, 1, by = 0.2),
                     limits = c(0, 1.1)) +
  guides(color = "none", linetype = "none") +
  theme_light()


fig4_6 <- simcirc2_DFPCA$fpc %>%
  mutate(ind = factor(1:n),
         gr = factor(rep(c(1, 2), each = n/2))) %>%
  ggplot(data = ., aes(x = X1, y = X2)) +
  geom_point(aes(color = gr, shape = gr)) +
  scale_color_manual(values=c("black", "tomato")) +
  xlab("PC1") + ylab("PC2") +
  guides(color = "none", shape = "none") +
  theme_light()


fig4_1
fig4_2
fig4_3
simcirc2_DFPCA$fig_mov1
simcirc2_DFPCA$fig_mov2
fig4_6

