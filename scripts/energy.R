# ==========================================================
# File: energy.R
# Purpose: Functions for analyzing the US electricity generation data (Fig 8)
# Author: Soobin Kim
# ==========================================================
source("R/compositional.R")

library(dplyr)
library(reshape2)
library(tidyr)
library(magrittr)
library(frechet)
library(fdapace)
library(ggplot2)
library(gridExtra)
library(ggrepel)
library(compositions)
library(Ternary)
library(cowplot)

# Load preprocessed data
load("data/energy.RData")

# Calculate ball statistics
energy_res <- list()
par(mfrow = c(3, 4), mar = c(2,3,2,3))
for(s in rownames(comp_y)){
  energy_res[[s]] <- 
    CompBMTraj(energy_dist, comp_y, x0 = s)
}
energy_res <- data.frame(do.call("rbind", energy_res))
colnames(energy_res) <- colnames(energy_comp_wide)[2:4]
energy_res$States <- rep(energy_comp_wide$STATE, each = 56)
energy_res$radius <- rep(seq(0, 1.1, length = 56), 50)


energy_res2 <- list()
for(i in 1:56){
  energy_res2[[i]] <- CompDistTraj(rad = seq(0, 1.1, length = 56)[i])
}
energy_d1 <- sapply(1:56, function(x){energy_res2[[x]]})
energy_d1_df <- data.frame(energy_d1) 
names(energy_d1_df) <- seq(0, 1.1, length = 56)
energy_d1_df$States <- energy_comp_wide$STATE

# FPCA of distance trajectory
e_DFPCA <- FPCA(
  Ly = data.frame(t(energy_d1)),
  Lt = rep(list(seq(0, 1.1, length = 56)), 50),
  optns = list(methodMuCovEst = "smooth", userBwCov = 0.05, userBwMu = 0.05)
)


# Plot
#### Figure 8
colors_50 <- createPalette(50, seed = c("#000000", "#FFFFFF"))
names(colors_50) <- unique(energy_res$States)

base_plot <- function() {
  par(mfrow = c(1,1), mar = c(0,0,0,0))
  TernaryPlot(
    alab = "Natural Gas", blab = "Other Fossil", clab = "Renewable / Nuclear",
    lab.offset = 0.12,
    tip.cex = 1, axis.cex = 0.7, tip.font = 1,
    grid.lines = 2,  grid.minor.lty = "dotted", grid.lty = "dotted",
    grid.col = "gray", 
    grid.lwd = 0.7, grid.minor.lwd = 0.7
  )
  AddToTernary(
    graphics::points,
    lapply(asplit(dplyr::filter(energy_res, radius == 0)[, 1:3], 1),
           function(x){as.numeric(x)}),
    col = colors_50,
    pch = 1, cex = 2, bg = "white")
  AddToTernary(
    text,
    lapply(asplit(dplyr::filter(energy_res, radius == 0)[, 1:3], 1),
           function(x){as.numeric(x)}),
    col = colors_50,
    unique(energy_res$States), cex = 0.6, font = 2)
  
}
# Convert base plot to a compatible object
fig8_1 <- ggdraw(base_plot)

fig8_2 <- reshape2::melt(energy_d1_df, id.vars = "States", value.name = "d", variable.name = "r") %>%
  mutate(r = as.numeric(levels(r))[r], d = d,
         label = if_else(r == min(r), as.character(States), NA_character_)) %>%
  ggplot(aes(x = r, y = d, group = States, colour = States)) +
  geom_line(alpha = 0.9) +
  scale_color_manual(values= colors_50,  
                     name = "States") +
  scale_x_continuous(breaks = seq(0, 1, by = 0.2),
                     limits = c(-0.01, 1.1)) +
  guides(color = "none") +
  ggrepel::geom_text_repel(aes(label = label),
                           size = 3, 
                           color = rep(colors_50, 56),
                           nudge_x = -0.01, direction = "y", hjust = "right",
                           na.rm = TRUE, max.overlaps = 15,
                           segment.size  = 0.1)+
  theme_light()

fig8_3 <- data.frame(e_DFPCA$xiEst[, 1:2], st = energy_d1_df$States) %>%
  ggplot(data = ., aes(x = X1, y = X2, label = st)) +
  geom_point(color = colors_50) +
  ggrepel::geom_text_repel(max.overlaps = 20, color = colors_50,
                           segment.size  = 0.1, size = 3) +
  xlab("FPC1") + ylab("FPC2") +
  theme_light()

fig8_4 <- data.frame(
  r = seq(0, 1.1, length = 56),
  mu = e_DFPCA$mu,
  mov1 = e_DFPCA$mu - sqrt(e_DFPCA$lambda[1]) * e_DFPCA$phi[, 1],
  mov1n = e_DFPCA$mu + sqrt(e_DFPCA$lambda[1]) * e_DFPCA$phi[, 1]
) %>%
  melt(data = ., id.vars = "r", value.name = "val", variable.name = "ftn") %>%
  ggplot(data = ., aes(x = r, y = val)) + 
  geom_line(aes(colour = ftn, linetype = ftn)) +
  scale_linetype_manual(
    name = "First mode (96.6%)",
    labels = c(
      expression("E[D"[X] * "(r)]"), 
      expression("E[D"[X] * "(r)]" * " - " * sqrt(lambda[1]) * phi[1](r)), 
      expression("E[D"[X] * "(r)]" * " + " * sqrt(lambda[1]) * phi[1](r))),
    values = c(1, 3, 5)) +
  scale_color_manual(
    name = "First mode (96.6%)",
    labels = c(
      expression("E[D"[X] * "(r)]"), 
      expression("E[D"[X] * "(r)]" * " - " * sqrt(lambda[1]) * phi[1](r)), 
      expression("E[D"[X] * "(r)]" * " + " * sqrt(lambda[1]) * phi[1](r))),
    values = c("black", "red", "blue")) +
  ylab("") +
  ylim(-0.01, 0.41) +
  theme_light() +
  theme(legend.position = "inside",
        legend.position.inside = c(0.77, 0.77))

fig8_5 <- data.frame(
  r = seq(0, 1.1, length = 56),
  mu = e_DFPCA$mu,
  mov1 =  e_DFPCA$mu - sqrt(e_DFPCA$lambda[2]) * e_DFPCA$phi[, 2],
  mov1n = e_DFPCA$mu + sqrt(e_DFPCA$lambda[2]) * e_DFPCA$phi[, 2]
) %>%
  melt(data = ., id.vars = "r", value.name = "val", variable.name = "ftn") %>%
  ggplot(data = ., aes(x = r, y = val)) + 
  geom_line(aes(colour = ftn, linetype = ftn)) +
  scale_linetype_manual(
    name = "Second mode (2.3%)",
    labels = c(
      expression("E[D"[X] * "(r)]"), 
      expression("E[D"[X] * "(r)]" * " - " * sqrt(lambda[1]) * phi[1](r)), 
      expression("E[D"[X] * "(r)]" * " + " * sqrt(lambda[1]) * phi[1](r))),
    values = c(1, 3, 5)) +
  scale_color_manual(
    name = "Second mode (2.3%)",
    labels = c(
      expression("E[D"[X] * "(r)]"), 
      expression("E[D"[X] * "(r)]" * " - " * sqrt(lambda[1]) * phi[1](r)), 
      expression("E[D"[X] * "(r)]" * " + " * sqrt(lambda[1]) * phi[1](r))),
    values = c("black", "red", "blue")) +
  ylab("") +
  ylim(-0.01, 0.41) +
  theme_light() +
  theme(legend.position = "inside",
        legend.position.inside = c(0.76, 0.77))

energy_mds <- reshape(energy_dist, idvar = "X", timevar = "Y", direction = "wide") %>%
  select(-X) %>%
  cmdscale()

fig8_6 <- data.frame(energy_mds) %>%
  mutate(States = energy_dist$X[1:50]) %>%
  ggplot(data = ., aes(x = X1, y = X2, label = States)) +
  geom_point(color = colors_50) +
  ggrepel::geom_text_repel(max.overlaps = 20, color = colors_50,
                           segment.size  = 0.1, size = 3) +
  xlab("MDS1") + ylab("MDS2") +
  theme_light()