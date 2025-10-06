# ==========================================================
# File: mortality.R
# Purpose: Functions for analyzing the human mortality data (Fig 6, 7, S2, S3)
# Author: Soobin Kim
# ==========================================================
source("R/distribution.R")

library(dplyr)
library(reshape2)
library(tidyr)
library(magrittr)
library(frechet)
library(fdapace)
library(ggplot2)
library(gridExtra)
library(ggrepel)

# Load preprocessed data
load("data/mortality.RData")

# Colors to distinguish EE vs rest
rb <- rainbow(41)
mort_colortbl <- data.frame(
  col = rb,
  country = country_names,
  cnt_full = c(
    "Australia", "Austria", "Belgium", "Bulgaria", "Belarus", "Canada",
    "Switzerland", "Chile", "Czechia", "Germany", "Denmark", "Spain",
    "Estonia", "Finland", "France", "U.K.", "Greece", "Hong Kong", "Croatia",
    "Hungary", "Ireland", "Iceland", "Israel", "Italy", "Japan", "Republic of Korea",
    "Lithuania", "Luxembourg", "Latvia", "Netherlands", "Norway",
    "New Zealand", "Poland", "Portugal", "Russia", "Slovakia", "Slovenia",
    "Sweden", "Taiwan", "Ukraine", "U.S.A."),
  ee = c(0, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0,
         0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 1, 0)
)
mort_colortbl$ee_br <- ifelse(mort_colortbl$ee == 1, "firebrick3", "darkgray")

# Calculate ball statistics
mortality2005 <- list()
for(cnt in country_names){
  mortality2005[[cnt]] <-
    DenBMTraj(ddf = mortality_dist, stck = mortality_stacked, x0 = cnt)
  print(cnt)
}

mortality2005_2 <- list()
for(i in 0:55){
  mortality2005_2[[i+1]] <- DenDistTraj(res = mortality2005, rad = i)
}
mortality_d1 <- sapply(1:56, function(x){mortality2005_2[[x]]})
mortality_d1_df <- data.frame(mortality_d1)
names(mortality_d1_df) <- seq(0, 1.1, length = 56)
mortality_d1_df$country <- country_names

# FPCA of distance trajectory
m_DFPCA <- FPCA(
  Ly = data.frame(t(mortality_d1)),
  Lt = rep(list(seq(0, 1.1, length = 56)), 41),
  optns = list(methodMuCovEst = "smooth", userBwCov = 0.05, userBwMu = 0.05)
)

# Transport: Separate pos, neg before integration
qSup <- mortality2005$AUS[[56]]$qSup
mort_T_pos <- mort_T_neg <- array(dim = c(41, 56))
for(i in 1:41){
  for(r in 1:56){
    tot <- (mortality2005[[i]][[r]]$qout - mortality2005[[i]][[56]]$qout)
    mort_T_pos[i,r] <- pracma::trapz(
      qSup[which(tot >= 0)], tot[which(tot >= 0)]
    )
    mort_T_neg[i,r] <- pracma::trapz(
      qSup[which(tot < 0)], tot[which(tot < 0)]
    )
  }
}

# FPCA of transport trajectory
mort_TFPCA_pos <- FPCA(
  Ly = data.frame(t(mort_T_pos)),
  Lt = rep(list(seq(0, 1.1, length = 56)), 41),
  optns = list(methodMuCovEst = "smooth", userBwCov = 0.05, userBwMu = 0.05)
)
mort_TFPCA_neg <- FPCA(
  Ly = data.frame(t(mort_T_neg)),
  Lt = rep(list(seq(0, 1.1, length = 56)), 41),
  optns = list(methodMuCovEst = "smooth", userBwCov = 0.05, userBwMu = 0.05)
)


# Plot
#### Figure 6
temp = replicate(3, data.frame(age = c(), density = c(), country = c()))
grid = seq(0, 1.1, length = 56)
use_rad = c(1, 26, 51)

for(j in 1:3){
  for(i in 1:41){
    temp[[j]] <- rbind(temp[[j]], data.frame(
      age = mortality2005[[i]][[use_rad[j]]]$dSup,
      density = mortality2005[[i]][[use_rad[j]]]$dout,
      country = rep(names(mortality2005)[i], 101)))
  }
  temp[[j]]$gr <- factor(ifelse(temp[[j]]$country %in% 
                     mort_colortbl$country[which(mort_colortbl$ee_br == "darkgray")],
                   "darkgray", "firebrick3"))
}

lapply(1:3, function(j){
  ggplot(data = temp[[j]], 
         aes(x = age, y = density, color = gr, group = country)) +
    geom_line(alpha = 0.7) +
    scale_color_manual(values=c("darkgray", "firebrick3")) +
    labs(title = paste("r =", grid[use_rad[j]])) +
    ylim(c(0, 0.044)) +
    guides(color = "none") +
    theme_light()
})


#### Figure 7
fig7_1 <- reshape2::melt(mortality_d1_df, id.vars = "country", value.name = "d", variable.name = "r") %>%
  mutate(r = as.numeric(levels(r))[r], d = d,
         label = if_else(r == min(r), as.character(country), NA_character_)) %>%
  mutate(label = if_else(label %in% mort_colortbl$country[which(mort_colortbl$ee ==1)], label, NA_character_)) %>%
  mutate(label = mort_colortbl$cnt_full[match(label, mort_colortbl$country)]) %>%
  ggplot(aes(x = r, y = d, group = country, colour = country)) +
  geom_line(alpha = 0.7) +
  scale_color_manual(
    values = mort_colortbl$ee_br, 
    name = "Country",
    labels = mort_colortbl$cnt_full) +
  scale_x_continuous(breaks = seq(0, 1, by = 0.2),
                     limits = c(-0.01, 1.1)) +
  guides(color = "none") +
  ggrepel::geom_text_repel(aes(label = label),
                           size = 3, color = "black",
                           nudge_x = -0.01, direction = "y", hjust = "right",
                           na.rm = TRUE, max.overlaps = 13,
                           segment.size  = 0.1) +
  theme_light()

fig7_2 <- data.frame(
  r = seq(0, 1.1, length = 56),
  mu = m_DFPCA$mu,
  mov1 = m_DFPCA$mu - sqrt(m_DFPCA$lambda[1]) * m_DFPCA$phi[, 1],
  mov1n = m_DFPCA$mu + sqrt(m_DFPCA$lambda[1]) * m_DFPCA$phi[, 1]
) %>%
  melt(data = ., id.vars = "r", value.name = "val", variable.name = "ftn") %>%
  ggplot(data = ., aes(x = r, y = val)) + 
  geom_line(aes(colour = ftn, linetype = ftn)) +
  scale_linetype_manual(
    name = "First mode (92.8%)",
    labels = c(
      expression("E[D"[X] * "(r)]"), 
      expression("E[D"[X] * "(r)]" * " - " * sqrt(lambda[1]) * phi[1](r)), 
      expression("E[D"[X] * "(r)]" * " + " * sqrt(lambda[1]) * phi[1](r))),
    values = c(1, 3, 5)) +
  scale_color_manual(
    name = "First mode (92.8%)",
    labels = c(
      expression("E[D"[X] * "(r)]"), 
      expression("E[D"[X] * "(r)]" * " - " * sqrt(lambda[1]) * phi[1](r)), 
      expression("E[D"[X] * "(r)]" * " + " * sqrt(lambda[1]) * phi[1](r))),
    values = c("black", "red", "blue")) +
  ylab("") +
  ylim(-0.01, 0.4) +
  theme_light() +
  theme(legend.position = "inside",
        legend.position.inside = c(0.77, 0.77))

fig7_3 <- data.frame(m_DFPCA$xiEst[, 1:2], country_names) %>%
  ggplot(data = ., aes(x = X1, y = X2, label = mort_colortbl$cnt_full)) +
  geom_point(aes(color = country_names, shape = country_names)) +
  scale_color_manual(
    values = mort_colortbl$ee_br,
    name = "Country", labels = mort_colortbl$cnt_full) +
  scale_shape_manual(
    values = ifelse(mort_colortbl$ee_br == "darkgray", 19, 17),
    name = "Country", labels = mort_colortbl$cnt_full) +
  ggrepel::geom_text_repel(max.overlaps = 30, color = "black",
                           segment.size  = 0.1) +
  xlab("PC1") + ylab("PC2") +
  guides(color = "none", shape = "none") +
  theme_light()

fig7_4 <- data.frame(
  r = seq(0, 1.1, length = 56),
  mu = m_DFPCA$mu,
  mov1 =  m_DFPCA$mu - sqrt(m_DFPCA$lambda[2]) * m_DFPCA$phi[, 2],
  mov1n = m_DFPCA$mu + sqrt(m_DFPCA$lambda[2]) * m_DFPCA$phi[, 2]
) %>%
  melt(data = ., id.vars = "r", value.name = "val", variable.name = "ftn") %>%
  ggplot(data = ., aes(x = r, y = val)) + 
  geom_line(aes(colour = ftn, linetype = ftn)) +
  scale_linetype_manual(
    name = "Second mode (4.2%)",
    labels = c(
      expression("E[D"[X] * "(r)]"), 
      expression("E[D"[X] * "(r)]" * " - " * sqrt(lambda[2]) * phi[2](r)), 
      expression("E[D"[X] * "(r)]" * " + " * sqrt(lambda[2]) * phi[2](r))),
    values = c(1, 3, 5)) +
  scale_color_manual(
    name = "Second mode (4.2%)",
    labels = c(
      expression("E[D"[X] * "(r)]"), 
      expression("E[D"[X] * "(r)]" * " - " * sqrt(lambda[2]) * phi[2](r)), 
      expression("E[D"[X] * "(r)]" * " + " * sqrt(lambda[2]) * phi[2](r))),
    values = c("black", "red", "blue")) +
  ylab("") +
  ylim(-0.01, 0.4) +
  theme_light() +
  theme(legend.position = "inside",
        legend.position.inside = c(0.76, 0.77))


#### Figure S2
MTPos <- as.data.frame(mort_T_pos) 
names(MTPos) <- seq(0, 1.1, length = 56)
MTPos$country <- country_names

figS2_1 <- melt(MTPos, id.vars = "country", value.name = "d", variable.name = "r") %>%
  mutate(r = as.numeric(levels(r))[r], d = d,
         label = if_else(r == min(r), as.character(country), NA_character_)) %>%
  mutate(label = if_else(label %in% mort_colortbl$country, label, NA_character_)) %>%
  mutate(label = mort_colortbl$cnt_full[match(label, mort_colortbl$country)]) %>%
  ggplot(aes(x = r, y = d, group = country, colour = country)) +
  ylim(0, 4) +
  geom_line(linewidth = 0.3, alpha = 0.7) +
  geom_abline(intercept = 0, slope = 0, color = "grey", linetype="dashed") +
  scale_color_manual(
    values=mort_colortbl$ee_br, name = "Country",
    labels = mort_colortbl$cnt_full) +
  scale_x_continuous(breaks = seq(0, 1, by = 0.2),
                     limits = c(-0.01, 1.1)) +
  labs(y = bquote({T[i]^"+"}(r))) + #, title = "Positive optimal transport trajectory") +
  ggrepel::geom_text_repel(aes(label = label),
                           size = 3, color = "gray50",
                           nudge_x = -0.01, direction = "y", hjust = "right",
                           na.rm = TRUE, max.overlaps = 5,
                           segment.size  = 0.1) +
  guides(color = "none") +
  theme_light()

figS2_2 <- data.frame(
  r = seq(0, 1.1, length = 56),
  mu = mort_TFPCA_pos$mu,
  mov1 = mort_TFPCA_pos$mu + sqrt(mort_TFPCA_pos$lambda[1]) * mort_TFPCA_pos$phi[, 1],
  mov1n = mort_TFPCA_pos$mu - sqrt(mort_TFPCA_pos$lambda[1]) * mort_TFPCA_pos$phi[, 1]
) %>%
  ggplot(data = ., aes(x = r)) +
  ylim(-0.05, 3) +
  geom_line(aes(y = mu)) +
  geom_line(aes(y = mov1), linetype = 1, col = "red") +
  geom_line(aes(y = mov1n), linetype = 1, col = "blue") +
  labs(title = "First mode of variation (94.2%)") + ylab("") +
  theme_light()

figS2_3 <- data.frame(mort_TFPCA_pos$xiEst[, 1:2], 
                      country_names = country_names[apply(mort_T_pos, 1, function(x){!all(round(x, 2) == 0)})]) %>%
  ggplot(data = ., aes(x = X1, y = -X2, label = mort_colortbl$cnt_full[apply(mort_T_pos, 1, function(x){!all(round(x, 2) == 0)})])) +
  geom_point(aes(color = country_names)) +
  scale_color_manual(
    values=mort_colortbl$ee_br[apply(mort_T_pos, 1, function(x){!all(round(x, 2) == 0)})],
    name = "Country",
    labels = mort_colortbl$cnt_full[apply(mort_T_pos, 1, function(x){!all(round(x, 2) == 0)})]
  ) +
  geom_text_repel(max.overlaps = 10, color = "gray50",
                  segment.size  = 0.1) +
  xlab("PC1") + ylab("PC2") +
  guides(color = "none") +
  theme(legend.spacing.y = unit(0.01, "lines")) +
  theme_light()

figS2_4 <- data.frame(
  r = seq(0, 1.1, length = 56),
  mu = mort_TFPCA_pos$mu,
  mov1 =  mort_TFPCA_pos$mu - sqrt(mort_TFPCA_pos$lambda[2]) * mort_TFPCA_pos$phi[, 2],
  mov1n = mort_TFPCA_pos$mu + sqrt(mort_TFPCA_pos$lambda[2]) * mort_TFPCA_pos$phi[, 2]
) %>%
  ggplot(data = ., aes(x = r)) +
  ylim(-0.05, 3) +
  geom_line(aes(y = mu)) +
  geom_line(aes(y = mov1), linetype = 1, col = "red") +
  geom_line(aes(y = mov1n), linetype = 1, col = "blue") +
  labs(title = "Second mode of variation (5.1%)") + ylab("") +
  theme_light()


grid.arrange(figS2_1, figS2_2, figS2_3, figS2_4, nrow = 2)

#### Figure S3
MTNeg <- as.data.frame(mort_T_neg) 
names(MTNeg) <- seq(0, 1.1, length = 56)
MTNeg$country <- country_names

figS3_1 <- melt(MTNeg, id.vars = "country", value.name = "d", variable.name = "r") %>%
  mutate(r = as.numeric(levels(r))[r], d = d,
         label = if_else(r == min(r), as.character(country), NA_character_)) %>%
  mutate(label = if_else(label %in% mort_colortbl$country, label, NA_character_)) %>%
  mutate(label = mort_colortbl$cnt_full[match(label, mort_colortbl$country)]) %>%
  ggplot(aes(x = r, y = d, group = country, colour = country)) +
  ylim(-9, 0) +
  geom_line(linewidth = 0.3, alpha = 0.7) +
  geom_abline(intercept = 0, slope = 0, color = "grey", linetype="dashed") +
  scale_color_manual(
    values=mort_colortbl$ee_br, name = "Country",
    labels = mort_colortbl$cnt_full) +
  scale_x_continuous(breaks = seq(0, 1, by = 0.2),
                     limits = c(-0.01, 1.1)) +
  labs(y = bquote({T[i]^"-"}(r))) + #, title = "Positive optimal transport trajectory") +
  ggrepel::geom_text_repel(aes(label = label),
                           size = 3, color = "gray50",
                           nudge_x = -0.01, direction = "y", hjust = "right",
                           na.rm = TRUE, max.overlaps = 5,
                           segment.size  = 0.1) +
  guides(color = "none") +
  theme_light()

figS3_2 <- data.frame(
  r = seq(0, 1.1, length = 56),
  mu = mort_TFPCA_neg$mu,
  mov1 =  mort_TFPCA_neg$mu - sqrt(mort_TFPCA_neg$lambda[1]) * mort_TFPCA_neg$phi[, 1],
  mov1n = mort_TFPCA_neg$mu + sqrt(mort_TFPCA_neg$lambda[1]) * mort_TFPCA_neg$phi[, 1]
) %>%
  ggplot(data = ., aes(x = r)) +
  ylim(-7.5, 0.5) +
  geom_line(aes(y = mu)) +
  geom_line(aes(y = mov1), linetype = 1, col = "red") +
  geom_line(aes(y = mov1n), linetype = 1, col = "blue") +
  labs(title = "First mode of variation (94.6%)") + ylab("") +
  theme_light()

figS3_3 <- data.frame(mort_TFPCA_neg$xiEst[, 1:2], 
                      country_names = country_names[mort_neg_col]) %>%
  ggplot(data = ., aes(x = -X1, y = X2, label = mort_colortbl$cnt_full[mort_neg_col])) +
  geom_point(aes(color = country_names)) +
  scale_color_manual(
    values=mort_colortbl$ee_br[mort_neg_col],
    name = "Country",
    labels = mort_colortbl$cnt_full[mort_neg_col]
  ) +
  geom_text_repel(max.overlaps = 10, color = "gray50",
                  segment.size  = 0.1) +
  xlab("PC1") + ylab("PC2") +
  guides(color = "none") +
  theme_light()

figS3_4 <- data.frame(
  r = seq(0, 1.1, length = 56),
  mu = mort_TFPCA_neg$mu,
  mov1 =  mort_TFPCA_neg$mu + sqrt(mort_TFPCA_neg$lambda[2]) * mort_TFPCA_neg$phi[, 2],
  mov1n = mort_TFPCA_neg$mu - sqrt(mort_TFPCA_neg$lambda[2]) * mort_TFPCA_neg$phi[, 2]
) %>%
  ggplot(data = ., aes(x = r)) +
  ylim(-7.5, 0.5) +
  geom_line(aes(y = mu)) +
  geom_line(aes(y = mov1), linetype = 1, col = "red") +
  geom_line(aes(y = mov1n), linetype = 1, col = "blue") +
  labs(title = "Second mode of variation (4.1%)") + ylab("") +
  theme_light()


grid.arrange(figS3_1, figS3_2, figS3_3, figS3_4, nrow = 2)
