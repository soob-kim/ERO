# ==========================================================
# File: mortality_preprocess.R
# Purpose: Functions for preprocessing the human mortality data
# Author: Soobin Kim
# ==========================================================
library(dplyr)
library(reshape2)
library(tidyr)
library(magrittr)
library(frechet)

# Load Data
file_names <- list.files("data/bltper_1x1")
mortality_raw <- list()
for(i in 1:length(file_names)){
  f <- paste("data/bltper_1x1/", file_names[i], sep = "")
  mortality_raw[[i]] <- read.table(f, skip = 2, header = T)
}
country_names <- sapply(
  strsplit(file_names, split = "[.]"),
  function(x){x[1]})
names(mortality_raw) <- country_names

# Remove regions
country_names <- setdiff(country_names, 
                         c("FRACNP", "DEUTE", "DEUTW", "NZL_MA", "NZL_NM",
                           "GBRTENW", "GBRCENW", "GBR_SCO", "GBR_NIR"))
mortality_raw <- mortality_raw[country_names]
names(mortality_raw) <- substr(country_names, 1, 3)
country_names <- substr(country_names, 1, 3)

# Easier-to-use form
mortality <- do.call(rbind.data.frame, mortality_raw)
mortality$Country <- rep(
  country_names, 
  sapply(country_names, 
         FUN = function(x){nrow(mortality_raw[[x]])}))
mortality$dx <- as.numeric(mortality$dx)
mortality$Age <- as.numeric(mortality$Age)
mortality$Age[is.na(mortality$Age)] <- 111
mortality <- select(mortality, Country, Year, Age, dx)


plot_density <- function(year, age_min = 40, age_max = 100, d = mortality_raw, plot = T){
  tt <- lapply(country_names, FUN = function(x){
    ttt <- d[[x]] %>%
      filter(Age != "110+") %>%
      mutate(Age = as.numeric(Age)) %>%
      filter(Year == year, Age <= age_max, Age >= age_min) %>%
      select(Age, dx)
    ttt$dx <- as.numeric(ttt$dx)
    if(nrow(ttt) == 0){
      return(NULL)
    }else{
      CreateDensity(freq = ttt$dx, bin = c(0:nrow(ttt)))
    }
  })
  names(tt) <- country_names
  tt <- tt[lengths(tt) != 0]
  if(plot){
    plot(tt[[1]]$x, tt[[1]]$y, type = "l", main = year,
         ylim = c(0, max(sapply(tt, function(x){max(x$y)}))),
         xlab = "Age", ylab = "Density")
    for(i in 2:length(tt)){
      lines(tt[[i]]$x, tt[[i]]$y)
    }
  }
  return(tt)
}

# Reshape data into frechet-usable form
stack_hist <- function(year, age_min = 40, age_max = 100, d = mortality_raw){
  tt <- lapply(country_names, FUN = function(x){
    ttt <- d[[x]] %>% 
      filter(Age != "110+") %>%
      mutate(Age = as.numeric(Age)) %>%
      filter(Year == year, Age <= age_max, Age >= age_min) %>%
      #rename(x = Age, y = dx) %>%
      select(Age, dx)
    ttt <- rep(ttt$Age, ttt$dx)
    if(length(ttt) == 0){
      return(NULL)
    }else{
      #hist(ttt, breaks = 90, plot = F)
      ttt
    }
    # if(nrow(ttt) == 0){
    #   return(NULL)
    # }
  })
  names(tt) <- country_names
  # tt <- tt[lengths(tt) != 0]
  # plot(tt[[1]]$x, tt[[1]]$y, type = "l", main = year,
  #      ylim = c(0, max(sapply(tt, function(x){max(x$y)}))),
  #      xlab = "Age", ylab = "Density")
  # for(i in 2:length(tt)){
  #   lines(tt[[i]]$x, tt[[i]]$y)
  # }
  tt <- tt[lengths(tt) != 0]
  return(tt)
}


# Get pairwise distance matrix (W2 distance)
get_distdf_den <- function(year, age_min, age_max, d = mortality_raw){
  use <- plot_density(year = year, age_min = age_min, age_max = age_max, plot = F)
  res = expand.grid(X = names(use), Y = names(use))
  nit = 1
  for(i in 1:length(use)){
    for(j in 1:length(use)){
      res$dist[nit] = dist4den(d1 = use[[i]], d2 = use[[j]])
      nit = nit + 1
    }
  }
  return(res)
}

# Save clean data
mortality_dist <- get_distdf_den(year = 2005, age_min = 40, age_max = 100)
mortality_stacked <- stack_hist(year = 2005)
save(country_names, mortality_raw, mortality, mortality_dist, mortality_stacked,
     file = "data/mortality.RData")
