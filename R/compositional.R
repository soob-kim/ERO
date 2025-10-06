# ==========================================================
# File: compositional.R
# Purpose: Functions for calculating ball statistics for Compositional data
# Author: Soobin Kim
# ==========================================================

library(frechet)
library(dplyr)
library(compositions)
library(ternary)

## Base functions: Most of them are copied & adjusted from frechet package
source(data/energy_preprocess.R)

#' Find the ball mean for Compositional data
#'
#' Args:
#'   distdf: A Dataframe of 3 columns (Obj 1, Obj 2, d(Obj 1, Obj 2))
#'   compobj: A matrix of size n * k (k-1 simplex)
#'   x0: Value of object X_i serving as the center
#'   radius: Numeric value specifying the radius of the ball.
#'     This should be a NOT scaled value.
#'
#' Returns: Ball mean (list)
CompBM <- function(distdf, compobj, x0, radius){
  tmp <- distdf %>%
    filter(X == x0) %>%
    filter(Z <= radius) 
  in_ball <- as.character(tmp$Y)
  mu_r <- CompFMean(yin = compobj[in_ball, , drop = F])
  return(mu_r)
}


#' Calculate ball mean trajectory \mu_i(r) over a grid of radii
#'
#' Computes the ball mean trajectory for a sequence of radii,
#' and the centers are each points in data
#'
#' Args:
#'   distdf: A Dataframe of 3 columns (Obj 1, Obj 2, d(Obj 1, Obj 2))
#'   compobj: A matrix of size n * k (k-1 simplex)
#'   x0: Value of object X_i serving as the center
#'   ngrid: Integer. Number of equally spaced grid points for the 
#'     radius at which the ball mean trajectory is calculated. 
#'     Default is 56.
#'   maxgrid: Numeric. Maximum radius for the grid. If not provided, it is 
#'     automatically set to 1.1.
#'   plot: Logical; Plot the results?
#'
#' Returns: A list of length ngrid; ith element is BM at ith radius at the grid

CompBMTraj <- function(distdf, compobj, x0, 
                       ngrid = 56, maxgrid = NULL,  plot = T){
  distdf_sc = distdf
  distdf_sc$Z = distdf_sc$Z / max(distdf_sc$Z)
  if(is.null(maxgrid)){
    maxgrid <- 1.1
  }
  res = matrix(nrow = ngrid, ncol = ncol(compobj))
  grid = seq(0, maxgrid, length = ngrid)
  
  for(i in 1:ngrid){
    res[i,] <- CompBM(distdf_sc, compobj, x0, diam = grid[i])$FMean
  }
  if(plot){
    colfunc <- colorRampPalette(c("red", "blue"))
    plot(compositions::acomp(res), col = colfunc(ngrid),
         pch = c(19, rep(1, ngrid-2), 19),
         labels = colnames(compobj))
    text(0.1, 0.5, x0)
    nlegend <- ifelse(ngrid>=10, 10, ngrid)
    textlegend <- grid[round(seq(from = ngrid, to = 1, length = nlegend))]
    textlegend <- round(textlegend, 2)
    plotrix::color.legend(
      1,0.2,1.2,1, textlegend,
      rev(colfunc(nlegend)), gradient = "y")
  }
  return(res)
}


#' Calculate distance trajectory D_i(r) 
#' 
#' Computes distance trajectory for a grid based on EucBMTraj output
#'
#' Args:
#'   res: Output (df) of `CompBMTraj`
#'   rad: Radius grid value
#'   ddf: Distance dataframe
#'
#' Returns: D: Vector of length n; distance trajectory
#' 
CompDistTraj <- function(res_df = energy_res, rad, ddf = energy_dist){
  k = dim(res_df) - 2
  #X_i(r)
  use <- dplyr::filter(res_df, radius == rad)
  #mu_oplus
  mu <- as.numeric(dplyr::filter(res_df, radius == 1)[1, 1:k])
  #d(X_i(r), mu_oplus)
  D <- apply(use[, 1:k], 1, function(x){CompGeoDist(x, mu)})
  if(scale){
    D <- D/max(ddf$Z)
  }
  return(D)
}



