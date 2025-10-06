# ==========================================================
# File: euclidean.R
# Purpose: Functions for calculating ball statistics for Euclidean data
# Author: Soobin Kim
# ==========================================================

library(dplyr)
library(fdapace)

#' Find the ball mean for Euclidean data
#'
#' Args:
#'   data: A numeric matrix of size n * p, where n is the number of 
#'     observations and p is the dimension.
#'   x0: Integer index of the row in `data` to be used as the center.
#'   radius: Numeric value specifying the radius of the ball.
#'     This should be a scaled value, within [0, 1].
#'
#' Returns: Ball mean, a vector of length p
EucBM <- function(data, d = NULL, x0 = NULL, radius){
  if(radius >= 1.5){
    warning("Radius could be too large; Use scaled value (diam = 1).")
  }
  maxd <- NULL
  # Calculate distance matrix d if only data is provided
  if(is.null(d)){
    d <- as.matrix(dist(data, diag = T))
    maxd <- max(d)
    d <- d / maxd
  }
  
  # Locate elements in the ball
  BallInd <- which(d[x0, ] <= radius)
  
  # Calculate F Mean
  mu_r <- apply(data[BallInd, ,drop = F], 2, mean)
  return(mu_r)
}

#' Calculate ball mean trajectory \mu_i(r) over a grid of radii
#'
#' Computes the ball mean trajectory for a sequence of radii,
#' and the centers are each points in data
#'
#' Args:
#'   data: A numeric matrix of size n * p, where n is the number of 
#'     observations and p is the dimension.
#'   ngrid: Integer. Number of equally spaced grid points for the 
#'     radius at which the ball mean trajectory is calculated. 
#'     Default is 56.
#'   maxr: Numeric. Maximum radius for the grid. If not provided, it is 
#'     automatically set to 1.1.
#'
#' Returns: A list of
#'   traj: (n \times ngrid) * 4 matrix; 
#'      Col 1:p : ball mean trajectory
#'      Col (p+1) : radius
#'      Col (p+2) : Index of center X_i
#'   maxd: Numeric; Maximum distance
#'   grid: Vector; grid of radii
EucBMTraj <- function(data, ngrid = 56, maxr = NULL){
  n <- dim(data)[1]
  p <- dim(data)[2]
  d <- as.matrix(dist(data, diag = T))
  maxd <- max(d)
  d <- d / maxd
  
  if(is.null(maxr)){
    maxr <- 1.1
  }
  # Setup grid
  grid <- seq(0, maxr, length = ngrid)
  
  res <- lapply(seq_len(n), function(x0) {
    mat <- matrix(NA, nrow = ngrid, ncol = p + 1)
    colnames(mat) <- c(paste(rep("X", p), 1:p, sep = ""), "r")
    for(i in 1:ngrid){
      mat[i, 1:p] <- EucBM(data = data, d = d, x0 = x0, radius = grid[i])
    }
    mat[, p + 1] <- grid   # fill last column
    mat
  })
  res_mat <- do.call(rbind, lapply(seq_along(res), function(i) {
    cbind(res[[i]], x0 = i)   # add x0 column
  }))
  
  return(list(traj = res_mat, maxd = maxd, grid = grid))
}

#' Calculate distance trajectory D_i(r) 
#' 
#' Computes distance trajectory for a grid based on EucBMTraj output
#'
#' Args:
#'   bmres: Output of `EucBMTraj`
#'
#' Returns: A list of
#'   traj: n * ngrid matrix; each column represents distance trajectory
#'   grid: Vector; grid of radii
EucDistTraj <- function(bmres){
  bmtraj <- bmres$traj
  maxd <- bmres$maxd
  grid <- bmres$grid
  
  n <- max(bmtraj[, "x0"])
  ngrid <- length(grid)
  p <- dim(bmtraj)[2] - 2
  
  res <- matrix(NA, nrow = n, ncol = ngrid)
  
  for(i in 1:ngrid){
    rad <- grid[i]
    
    # Calculate global Frechet mean mu_oplus
    mu <- apply(bmtraj[which(bmtraj[, "r"] == 0), 1:p], 2, mean)
    # Select \mu_i(rad) for given radius rad
    mu_i <- bmtraj[which(bmtraj[, "r"] == rad), 1:p]
    
    # Calculate D_i(r) = d(X_i(r), mu_oplus)
    D <- apply(mu_i, 1, function(x){norm(x - mu, type = "2")})
    D <- D / maxd # scale 
    res[, i] <- D
  }
  
  return(list(traj = res, grid = grid))
}

#' Wrapper function for FPCA of distance trajectories
#' 
#' Extracts information for plotting modes of variation and FPC scatterplot
#'
#' Args:
#'   dtres: Output of `EucDistTraj`
#'   ind: Index (names) of data observations
#'   opts: A list of options control parameters for fdapace::FPCA
#'   plot: Logical; produce modes of variation (MOV) plots?
#'   legend_loc: Vector of length 2; where to put legends in MOV plots?
#'   yliml: Numeric; lower limit for ylim in MOV plots
#'   ylimu: Numeric; upper limit for ylim in MOV plots
#'
#' Returns: A list of
#'  fve: Percentage of variance explained by the first two FPCs
#'  mov1: ngrid * 4 DF; r (radius grid), mu (mean function), 
#'      mov (1st MOV with alpha = -1), movn (1st MOV with alpha = 1)
#'  mov1: ngrid * 4 DF; r (radius grid), mu (mean function), 
#'      mov (2nd MOV with alpha = -1), movn (2nd MOV with alpha = 1)
#'  fpc: n * 2 DF; First two FPC scores
#'  fig_mov1: ggplot object for 1st MOV
#'  fig_mov2: ggplot object for 2nd MOV
EucDFPCA <- function(dtres, ind, 
                     opts = list(methodMuCovEst = "smooth", userBwCov = 0.05, userBwMu = 0.05),
                     plot = T, 
                     legend_loc = c(0.77, 0.77),
                     yliml = -0.01, ylimu = 0.35){
  n <- length(ind)
  
  # FPCA
  DFPCA <- FPCA(
    Ly = data.frame(t(dtres$traj)),
    Lt = rep(list(dtres$grid), n),
    optns = opts
  )
  
  # % of variance explained
  fve <- c(DFPCA$cumFVE[1],
           DFPCA$cumFVE[2] - DFPCA$cumFVE[1]) * 100
  
  # FPC scores
  fpc <- data.frame(DFPCA$xiEst[, 1:2], ind = ind)
  
  # Modes of variation
  mov1 <- data.frame(
    r = dtres$grid, mu = DFPCA$mu,
    mov =  DFPCA$mu - sqrt(DFPCA$lambda[1]) * DFPCA$phi[, 1],
    movn = DFPCA$mu + sqrt(DFPCA$lambda[1]) * DFPCA$phi[, 1]
  )
  
  mov2 <- data.frame(
    r = dtres$grid, mu = DFPCA$mu,
    mov =  DFPCA$mu - sqrt(DFPCA$lambda[2]) * DFPCA$phi[, 2],
    movn = DFPCA$mu + sqrt(DFPCA$lambda[2]) * DFPCA$phi[, 2]
  )
  
  # MOV plots
  fig_mov1 <- NULL
  fig_mov2 <- NULL
  if(plot){
    fig_mov1 <- mov1 %>%
      melt(data = ., id.vars = "r", value.name = "val", variable.name = "ftn") %>%
      ggplot(data = ., aes(x = r, y = val)) + 
      geom_line(aes(colour = ftn, linetype = ftn)) +
      scale_linetype_manual(
        name = paste("First mode (", round(fve, 1)[1], "%)", sep = ""),
        labels = c(
          expression("E[D"[X] * "(r)]"), 
          expression("E[D"[X] * "(r)]" * " - " * sqrt(lambda[1]) * phi[1](r)), 
          expression("E[D"[X] * "(r)]" * " + " * sqrt(lambda[1]) * phi[1](r))),
        values = c(1, 3, 5)) +
      scale_color_manual(
        name = paste("First mode (", round(fve, 1)[1], "%)", sep = ""),
        labels = c(
          expression("E[D"[X] * "(r)]"), 
          expression("E[D"[X] * "(r)]" * " - " * sqrt(lambda[1]) * phi[1](r)), 
          expression("E[D"[X] * "(r)]" * " + " * sqrt(lambda[1]) * phi[1](r))),
        values = c("black", "red", "blue")) +
      ylab("") + ylim(yliml, ylimu) +
      theme_light() +
      theme(legend.position = "inside",
            legend.position.inside = legend_loc)
    
    fig_mov2 <- mov2 %>%
      melt(data = ., id.vars = "r", value.name = "val", variable.name = "ftn") %>%
      ggplot(data = ., aes(x = r, y = val)) + 
      geom_line(aes(colour = ftn, linetype = ftn)) +
      scale_linetype_manual(
        name = paste("Second mode (", round(fve, 1)[2], "%)", sep = ""),
        labels = c(
          expression("E[D"[X] * "(r)]"), 
          expression("E[D"[X] * "(r)]" * " - " * sqrt(lambda[1]) * phi[1](r)), 
          expression("E[D"[X] * "(r)]" * " + " * sqrt(lambda[1]) * phi[1](r))),
        values = c(1, 3, 5)) +
      scale_color_manual(
        name = paste("Second mode (", round(fve, 1)[2], "%)", sep = ""),
        labels = c(
          expression("E[D"[X] * "(r)]"), 
          expression("E[D"[X] * "(r)]" * " - " * sqrt(lambda[1]) * phi[1](r)), 
          expression("E[D"[X] * "(r)]" * " + " * sqrt(lambda[1]) * phi[1](r))),
        values = c("black", "red", "blue")) +
      ylab("") + ylim(yliml, ylimu) + 
      theme_light() +
      theme(legend.position = "inside",
            legend.position.inside = legend_loc)
  }
  return(list(fve = fve, mov1 = mov1, mov2 = mov2, fpc = fpc,
              fig_mov1 = fig_mov1, fig_mov2 = fig_mov2))
}

