#'
#' @title Run the turnover approaches on the Jena data
#' 
#' @description Uses the AIC and SES - based turnover approaches to 
#' calculate the proportion of the species pool required to support different
#' ecosystem functions on the Jena data
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

# get the species and function names for the different datasets

# Jena data

# get the species names
jena_cols <- names(jena_dat)
jena_sp <- jena_cols[grepl(pattern = "[A-Z][a-z]{2}[.][a-z]{3}", jena_cols)]

# get the function names
jena_funcs <- c("biomass", "plantCN", "soilC", "soilorgC", "herbi", "micBMC", "phoact", "poll","rootBM")

# convert abundance data into presence absence
jena_dat[,jena_sp] <- vegan::decostand(jena_dat[,jena_sp], method = "pa")


# set the number of randomisations to do for the null expectation
n_null <- 999

# get the observed proportion of the species pool
aic_obs <- 
  prop_species_pool(data = jena_dat, 
                    func_names = jena_funcs, 
                    sp_names = jena_sp, method = "AIC")

# get the null proportion of the species pool
aic_ran <- 
  prop_species_pool_random(data = jena_dat, 
                           func_names = jena_funcs, 
                           sp_names = jena_sp, method = "AIC", n = n_null)

# combine into a list
aic_dat <- list( list(obs = aic_obs, null = aic_ran) )

# save as a RDS file
saveRDS(object = aic_dat, file = "code/02-paper-1/02-turnover-approaches/AIC_output.rds")

# run the SES-based turnover approach

# get the observed proportion of the species pool
ses_obs <- 
  prop_species_pool(data = jena_dat, 
                    func_names = jena_funcs, 
                    sp_names = jena_sp, method = "SES", n_ran = 999)

# get the null proportion of the species pool
ses_ran <- 
  prop_species_pool_random(data = jena_dat, 
                           func_names = jena_funcs, 
                           sp_names = jena_sp, method = "SES", n_ran = 999, 
                           n = n_null)

# combine into a list
ses_dat <- list( list(obs = ses_obs, null = ses_ran) )

# save as an RDS file
saveRDS(object = aic_dat, file = "code/02-paper-1/02-turnover-approaches/SES_output.rds")

### END
