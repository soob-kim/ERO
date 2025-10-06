# ==========================================================
# File: simCircOut.R
# Purpose: Simulation in circular data (Supplement Fig 1)
# Author: Soobin Kim
# ==========================================================
# - Sample size: n = 100
# - Distribution:
#     96 of theta = seq(pi/4, 3pi/4), seq(5pi/4, 7pi/4) 
#     4 of theta_out = 0, 0, pi, pi
#     rho \sim N(2, 0.2^2)

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
theta = c(seq(pi/4, 3*pi/4, length = 48),
          seq(5*pi/4, 7*pi/4, length = 48),
          0, 0, pi, pi)
rho = rnorm(n, mean = 1, sd = 0.1)

simcirc3 = cbind(2 * rho * cos(theta), 2 * rho * sin(theta))

# Calculate ball statistics
simcirc3_BMT <- EucBMTraj(simcirc3)
simcirc3_DT <- EucDistTraj(simcirc3_BMT)

# FPCA of distance trajectory
simcirc3_DFPCA <- EucDFPCA(simcirc3_DT, ind = factor(1:n),
                          legend_loc = c(0.23, 0.23),
                          yliml = -0.01, ylimu = 0.5)

# Plot
#### Figure 3
figS1_1 <- data.frame(simcirc3, gr = factor(c(rep(1, 96), rep(2, 4)))) %>%
  ggplot(data = ., aes(x = X1, y = X2, color = gr)) +
  geom_point() +
  scale_color_manual(values=c("black", "tomato")) +
  guides(color = "none") +
  theme_light()

figS1_2 <- plot_ly(type = 'scatter3d') 
figS1_2 <- add_trace(
  figS1_2, mode = "lines",
  x = simcirc3_BMT$traj[, "r"], 
  y = simcirc3_BMT$traj[, "X1"], 
  z = simcirc3_BMT$traj[, "X2"],
  color = as.factor(simcirc3_BMT$traj[, "x0"]),
  line = list(color = "gray", opacity = 1, width = 0.7), showlegend = F)
figS1_2 <- add_trace(
  figS1_2, mode = "markers", x = rep(0, n), 
  y = simcirc3[,1],
  z = simcirc3[,2],
  color = as.factor(simcirc3_BMT$traj[which(simcirc3_BMT$traj[, "r"] == 0), "x0"]),
  marker = list(size = 3, color = c(rep("midnightblue", 96), rep("tomato", 4))), 
  showlegend = F) %>%
  layout(scene = list(xaxis = list(title = 'r'),
                      yaxis = list(title = 'X1'),
                      zaxis = list(title = 'X2'),
                      aspectratio = list(x=1, y=1, z=1)))



figS1_3 <- as_tibble(simcirc3_DT$traj, 
                    .name_repair = ~ as.character(simcirc3_DT$grid)) %>%
  mutate(ind = 1:n) %>%
  melt(id.vars = "ind", value.name = "d", variable.name = "r") %>%
  mutate(r = as.numeric(levels(r))[r], d = d,
         ind = as.numeric(ind),
         gr = factor(if_else(ind <= 96, 1, 2))) %>%
  ggplot(aes(x = r, y = d, group = ind, colour = gr)) +
  ylim(0, 0.6) +
  geom_line(alpha = 0.5) +
  scale_color_manual(values=c("black", "tomato")) +
  scale_x_continuous(breaks = seq(0, 1, by = 0.2),
                     limits = c(0, 1.1)) +
  guides(color = "none") +
  theme_light()


figS1_6 <- simcirc3_DFPCA$fpc %>%
  mutate(gr = factor(c(rep(1, 96), rep(2, 4)))) %>%
  ggplot(data = ., aes(x = X1, y = X2)) +
  geom_point(aes(color = gr)) +
  scale_color_manual(values=c("black", "tomato")) +
  xlab("PC1") + ylab("PC2") +
  guides(color = "none") +
  theme_light()


figS1_1
figS1_2
figS1_3
simcirc3_DFPCA$fig_mov1
simcirc3_DFPCA$fig_mov2
figS1_6
