# ==========================================================
# File: energy_preprocess.R
# Purpose: Functions for preprocessing the US electricity generation
# Author: Soobin Kim
# ==========================================================
#### Load data
library(readxl)
energy <- read_xls("data/annual_generation_state.xls", col_names = F)
colnames(energy) <- energy[2,]

energy <- energy[-c(1,2), ]

#### Preprocess
library(dplyr)
energy_comp <- energy %>%
  mutate(GENERATION = as.numeric(`GENERATION (Megawatthours)`)) %>%
  select(YEAR, STATE, `ENERGY SOURCE`, GENERATION) %>%
  filter(YEAR == "2000", STATE != "US-TOTAL") %>%
  filter(!`ENERGY SOURCE` %in% c("Total", "Pumped Storage")) %>%
  mutate(ES_GR = ifelse(
    `ENERGY SOURCE` == "Natural Gas", "NatGas",
    ifelse(
      `ENERGY SOURCE` %in% c("Coal", "Petroleum", "Other Gases"), "OthFoss",
      "RenNuc"
    )
  )) %>%
  group_by(STATE, ES_GR) %>%
  summarise(GEN_SUM = sum(GENERATION)) %>%
  mutate(GEN_P = GEN_SUM / sum(GEN_SUM)) %>%
  select(STATE, ES_GR, GEN_P)
# excluded the “pumped storage” category 
# Natural Gas: “natural gas” alone
# Other Fossil: “coal”, “petroleum”, “other gases”
# Renewables and Nuclear: “hydroelectric conventional”, “solar thermal and photovoltaic”, “geothermal”, “wind”, “wood and wood derived fuels”, “other biomass”, “nuclear” and “other”

energy_comp <- filter(energy_comp, STATE != "DC")

library(reshape2)
energy_comp_wide <-
  dcast(energy_comp, STATE ~ ES_GR, value = c("GEN_P"))

energy_comp_wide[is.na(energy_comp_wide)] <- 0



# L2 norm
l2norm <- function(x){
  #sqrt(sum(x^2))
  as.numeric(sqrt(crossprod(x)))
}

CompGeoDist <- function(y1, y2) {
  if (abs(length(y1) - length(y2)) > 0) {
    stop("y1 and y2 should be of the same length.")
  }
  if ( !all(round(sum(y1), 5) == 1) ) {
    stop("y1 is not compositional data.")
  }
  if ( !all(round(sum(y2), 5) == 1) ) {
    stop("y2 is not compositional data.")
  }
  y1 = y1 / sum(y1)
  y2 = y2 / sum(y2)
  sqrt_y1 = sqrt(y1)
  sqrt_y2 = sqrt(y2)
  if (all(y1 == y2)){
    return(0)
  } else if (sum(sqrt_y1 * sqrt_y2) > 1){
    return(0)
  } else if (sum(sqrt_y1 * sqrt_y2) < -1){
    return(pi)
  } else return(acos(sum(sqrt_y1 * sqrt_y2)))
}

SpheGeoHess <- function(x,y) { #,tol = 1e-10){
  return(- sum(x * y) * (1 - sum(x * y) ^ 2) ^ (-1.5) * x %*% t(x))
}


SpheGeoGrad <- function(x,y) { 
  tmp <- 1 - sum(x * y) ^ 2
  return(- (tmp) ^ (-0.5) * x)
  # if (tmp < tol) {
  #   return(- Inf * x)
  # } else {
  #   return(- (tmp) ^ (-0.5) * x)
  # }
}
GloSpheGeoReg <- function(xin, yin, xout) {
  k = length(xout)
  n = length(xin)
  m = ncol(yin)
  
  xbar <- colMeans(xin)
  Sigma <- cov(xin) * (n-1) / n
  if (sum(abs(xin - 1)) == 0) {
    # compute the Fréchet mean
    invSigma <- diag(1, nrow(Sigma))
  } else {
    invSigma <- solve(Sigma)
  }
  
  yout = sapply(1:k, function(j){
    s <- 1 + t(t(xin) - xbar) %*% invSigma %*% (xout[j,] - xbar)
    s <- as.vector(s)
    
    # initial guess
    y0 = colMeans(yin*s)
    y0 = y0 / l2norm(y0)
    if ( any( sapply( 1:n, function(i) isTRUE( all.equal( sum(yin[i,]*y0), 1 ) ) ) ) ){
      # if (sum(sapply(1:n, function(i) sum(yin[i,]*y0)) > 1-1e-8)){
      #if (sum( is.infinite (sapply(1:n, function(i) (1 - sum(yin[i,]*y0)^2)^(-0.5) )[ker((xout[j] - xin) / bw)>0] ) ) + 
      #   sum(sapply(1:n, function(i) 1 - sum(yin[i,] * y0)^2 < 0)) > 0){
      y0[1] = y0[1] + 1e-3
      y0 = y0 / l2norm(y0)
    }
    
    objFctn = function(y){
      if ( ! isTRUE( all.equal(l2norm(y),1) ) ) {
        return(list(value = Inf))
      }
      f = mean(sapply(1:n, function(i) SpheGeoDist(yin[i,], y)^2))
      #f = mean(s * sapply(1:n, function(i) SpheGeoDist(yin[i,], y)^2))
      g = 2 * colMeans(t(sapply(1:n, function(i) SpheGeoDist(yin[i,], y) * SpheGeoGrad(yin[i,], y))) * s)
      res = sapply(1:n, function(i){
        grad_i = SpheGeoGrad(yin[i,], y)
        return((grad_i %*% t(grad_i) + SpheGeoDist(yin[i,], y) * SpheGeoHess(yin[i,], y)))
        #return((grad_i %*% t(grad_i) + SpheGeoDist(yin[i,], y) * SpheGeoHess(yin[i,], y)) * s[i])
      }, simplify = "array")
      h = 2 * apply(res, 1:2, mean)
      return(list(value=f, gradient=g, hessian=h))
    }
    res = trust::trust(objFctn, y0, 0.1, 1)
    return(res$argument)
  })
  return(t(yout))
}

CompFMean <- function(yin=NULL){
  if (!is.matrix(yin) | !is.numeric(yin))
    stop("yin should be a numerical matrix.")
  if ( !all(round(rowSums(yin), 5) == 1) ){
    yin = yin / rowSums(yin^2)
    warning("Each row of yin has been standardized to enforce sum equal to 1.")
  }
  n = nrow(yin)
  yout <- GloSpheGeoReg(xin = matrix(1, nrow=n), yin = sqrt(yin), xout = matrix(1, nrow=n))[1,]
  yout = yout^2
  res <- list(FMean = yout, FMean_acomp = compositions::acomp(yout))
  return(res)
}

get_distdf <- function(compobj, metric = CompGeoDist, id){
  ## compobj in matrix type for now
  res = expand.grid(X = id, Y = id)
  n = nrow(compobj)
  nit = 1
  for(i in 1:n){
    for(j in 1:n){
      res$Z[nit] = metric(compobj[i,], compobj[j,])
      nit = nit + 1
    }
  }
  return(res)
}

comp_y <- as.matrix(energy_comp_wide[, 2:4])
rownames(comp_y) <- energy_comp_wide$STATE

energy_dist <- get_distdf(
  comp_y, metric = CompGeoDist, id = energy_comp_wide$STATE)

save(energy, energy_comp, energy_comp_wide, comp_y, energy_dist,
     file = "data/energy.RData")