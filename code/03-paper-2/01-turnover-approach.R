#'
#' @title Run the turnover approaches on the Jena and BIODEPTH data
#' 
#' @description Uses the AIC and SES - based turnover approaches to 
#' calculate the proportion of the species pool required to support different
#' ecosystem functions on the Jena data and two BIODEPTH datasets: Sweden
#' and Portugal
#'

# load relevant libraries
library(dplyr)
library(ggplot2)
library(ggbeeswarm)
library(cowplot)

# load relevant scripts
source("code/helper-plotting-theme.R")
source("code/helper-turnover-approach-functions.R")

# load the relevant datasets
jena_dat <- readRDS("data/jena_data.rds")
swe_dat <- readRDS("data/biodepth_swe_data.rds")
prt_dat <- readRDS("data/biodepth_prt_data.rds")

# get the species and function names for the different datasets

# Jena data

# get the species names
jena_cols <- names(jena_dat)
jena_sp <- jena_cols[grepl(pattern = "[A-Z][a-z]{2}[.][a-z]{3}", jena_cols)]

# get the function names
jena_funcs <- c("biomass", "plantCN", "soilC", "soilorgC", "herbi", "micBMC", "phoact", "poll","rootBM")

# BIODEPTH: Sweden
swe_cols <- names(swe_dat)
swe_sp <- swe_cols[grepl(pattern = "[A-Z]{6}[1]", swe_cols)]

# get the function names
swe_funcs <- c("biomassY3", "root3", "cotton3", "wood3", "N.g.m2", "light3")

# BIODEPTH: Portugal
prt_cols <- names(prt_dat)
prt_sp <- prt_cols[grepl(pattern = "[A-Z]{6}[1]", prt_cols)]

# get the function names
prt_funcs <- c("biomassY3", "root3", "cotton3", "wood3", "N.Soil", "N.g.m2")

# pull these datasets into a list
dat_list <- list(jena_dat, swe_dat, prt_dat)
dat_funcs <- list(jena_funcs, swe_funcs, prt_funcs)
dat_sp <- list(jena_sp, swe_sp, prt_sp)

# convert abundance data into presence absence
dat_list <- 
  
  mapply(function(D, S){
  
  D[,S] <- vegan::decostand(D[,S], method = "pa")
  return(D)
  
}, dat_list, dat_sp)


# run the AIC turnover approach on these datasets
aic_dat <- 
  
  mapply(function(D, Func, S) {
  
  # get the observed proportion of the species pool
  aic_obs <- 
    prop_species_pool(data = D, 
                      func_names = Func, 
                      sp_names = S, method = "AIC")
  
  # get the null proportion of the species pool
  aic_ran <- 
    prop_species_pool_random(data = D,
                             func_names = Func,
                             sp_names = S, method = "AIC", n = 5)
  
  return( list(aic_obs, aic_ran) )
  
}, dat_list, dat_funcs, dat_sp, SIMPLIFY = FALSE)

# set set the names of the datasets
names(aic_dat) <- c("jena", "swe", "prt")

# run the SES-based turnover approach
ses_dat <- 
  
  mapply(function(D, Func, S) {
    
    # get the observed proportion of the species pool
    ses_obs <- 
      prop_species_pool(data = D, 
                        func_names = Func, 
                        sp_names = S, method = "SES", n_ran = 999)
    
    # get the null proportion of the species pool
    ses_ran <- 
      prop_species_pool_random(data = D,
                               func_names = Func,
                               sp_names = S, method = "SES", n_ran = 999, 
                               n = 5)
    
    return( list(ses_obs, ses_ran) )
    
  }, dat_list, dat_funcs, dat_sp, SIMPLIFY = FALSE)

# set set the names of the datasets
names(ses_dat) <- c("jena", "swe", "prt")



x <- SES_sp(data = jena_dat, func_names = jena_funcs, sp_names = jena_sp, n_ran = 100)








