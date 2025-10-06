# ==========================================================
# File: distribution.R
# Purpose: Functions for calculating ball statistics for Distributional data
# Author: Soobin Kim
# ==========================================================

library(frechet)
library(dplyr)

#' Find the ball mean for Distributional data
#'
#' Args:
#'   distdf: A Dataframe of 3 columns (Obj 1, Obj 2, d(Obj 1, Obj 2))
#'   stackobj: A list of distributional data (`stack_hist()`)
#'   x0: Value of object X_i serving as the center
#'   radius: Numeric value specifying the radius of the ball.
#'     This should be a NOT scaled value.
#'
#' Returns: Ball mean, frechet::DenFMean object
DenBM <- function(distdf, stackobj, x0, radius){
  tt <- distdf %>%
    filter(X == x0) %>%
    filter(dist <= radius) 
  in_ball <- as.character(tt$Y)
  mu_r <- DenFMean(yin = stackobj[in_ball])
  return(mu_r)
}



plot_BallDen <- function(res, ngrid = 56, x0, year, ylim = c(0, 0.045),
                         title = NULL){
  if(is.null(title)){
    title <- paste(x0, year)
  }
  ymax = sapply(res, function(x){max(x$dout)})
  colfunc <- colorRampPalette(c("gold", "red3"))
  cols <- colfunc(ngrid)
  plot(res[[1]]$dSup, res[[1]]$dout, type = "l", ylim = ylim, 
       col = cols[1], xlab = "Age", ylab = "Density", main = title)
  for(i in 2:ngrid){
    lines(res[[i]]$dSup, res[[i]]$dout, col = cols[i])
  }
  nlegend <- ifelse(ngrid>=12, 12, ngrid)
  grid = seq(0, 1.1, length = ngrid)
  textlegend <- grid[round(seq(from = ngrid, to = 1, length = nlegend))]
  textlegend <- round(textlegend, 2)
  # plotrix::color.legend(45,ylim[2]/2, 50,ylim[2], textlegend, align = "rb",
  #                       rev(colfunc(nlegend)), gradient = "y")
  PlotTools::SpectrumLegend(
    #"topright",
    x = 45, y = ylim[2],
    cex = 0.8, lwd = 15,
    palette = rev(colfunc(ngrid)),
    legend = seq(0, 1, length = 5),
    bty = "n", # No framing box
    xpd = NA, # Don't clip at margins
    # title.font = 2, # Bold. Supported from R 3.6 onwards
    title = "Radius"
  )
}

#' Calculate ball mean trajectory \mu_i(r) over a grid of radii
#'
#' Computes the ball mean trajectory for a sequence of radii,
#' and the centers are each points in data
#'
#' Args:
#'   ddf: A Dataframe of 3 columns (Obj 1, Obj 2, d(Obj 1, Obj 2))
#'   stck: A list of distributional data (`stack_hist()`)
#'   x0: Value of object X_i serving as the center
#'   ngrid: Integer. Number of equally spaced grid points for the 
#'     radius at which the ball mean trajectory is calculated. 
#'     Default is 56.
#'   maxgrid: Numeric. Maximum radius for the grid. If not provided, it is 
#'     automatically set to 1.1.
#'   plot: Logical; Plot the results?
#'
#' Returns: A list of length ngrid; ith element is BM at ith radius at the grid
DenBMTraj <- function(ddf, stck, x0, ngrid = 56, maxgrid = NULL, plot = T){
  res <- list()
  if(is.null(maxgrid)){
    maxgrid <- max(ddf$dist)*1.1
  }
  grid = seq(0, maxgrid, length = ngrid)
  for(i in 1:ngrid){
    res[[i]] <- DenBM(distdf = ddf, stackobj = stck, 
                      x0 = x0, radius = grid[i])
  }
  if(plot){
    plot_BallDen(res=res, ngrid=ngrid, x0=x0, year=year)
  }
  return(res)
}


#' Calculate distance trajectory D_i(r) 
#' 
#' Computes distance trajectory for a grid based on EucBMTraj output
#'
#' Args:
#'   res: Output of `DenBMTraj`
#'   rad: Index of radii grid
#'   ddf: Distance dataframe
#'
#' Returns: D: Vector of length n; distance trajectory
DenDistTraj <- function(res = mortality2005, rad, ddf = mortality_dist){
  ngrid = length(res)
  # X_i(r)
  use <- lapply(res, function(x){x[[rad + 1]]}) #length = country
  
  #d(X_i(r), mu_oplus)
  mu <- res[[1]][[ngrid]]
  D <- sapply(use, function(obj){
    dist4den(
      d1 = list(x = obj$dSup, y = obj$dout),
      d2 = list(x = mu$dSup, y = mu$dout))
  })
  D <- D/max(ddf$dist)
  return(D)
}

