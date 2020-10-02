
# Project: Mutifunctionality workshop (Lund 2019)

# Title: Preparing the model data to add multifunctionality into Lotka-Volterra models

# load libraries
library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)
library(here)
library(corrplot)

# choose scripts to draw functions from
source(here("Scripts/MF_functions_collated.R"))
source(here("Scripts/function_plotting_theme.R"))
source(here("Scripts/Multifunctionality-Simulations/Multifunc_simulations_functions.R"))

# load the Lotka-Volterra simulated data
lv_dat <- 
  read_csv(file = here("data/stachova_leps_model_mf.csv"))

# check dimensions of the data
nrow(lv_dat)

# make a vector of species names
sp_names <- unique(lv_dat$species)


# simulate three sets of species-specific functions

lv_mf <- 
  lapply(split(lv_dat, lv_dat$run), function(x) {
  
  # for each run that is applied to, generate functions with three correlation structures  
    
  # set correlation cluster values for each run
  lc <- 
    list(c1 = c(0, 0, 0),
         c2 = c(-0.4, -0.3, -0.5),
         c3 = c(0.1, 0.4, 0.4))
  
  # set an output list to gather the species-specific function values
  ss_funcs <- vector("list", length = length(lc))
  
  for (i in 1:length(lc)) {
    
    # choose the number of functions to simulate (cannot be varied without varying the correlation matrix)
    func_n <- 9
    
    # set the function names
    func_names <- paste0("F_", 1:func_n)
    
    # create an empty data matrix
    Funcs <- 
      matrix(nrow = length(sp_names),
             ncol = func_n,
             dimnames = list(sp_names, func_names))
    
    # number of species
    specnum <- length(sp_names)
    
    # choose pairwise correlation strength
    COR <- 0
    
    # make correlation matrix (strictly speaking a covariance matrix but for these simulations it does not matter)
    Sigma <- matrix(COR, ncol = func_n, nrow = func_n)
    
    # make three 'cluster' of correlated functions
    Sigma[1:3,1:3] <- lc[[i]][1]
    Sigma[8:9,8:9] <- lc[[i]][2]
    Sigma[6:7,6:7] <- lc[[i]][3]
    
    diag(Sigma) <- 1
    
    is.positive.definite(x = Sigma)
    
    # draw correlated functions (with mean 0)
    corF <- mvrnorm(n = specnum, mu = rep(0, func_n), Sigma = Sigma)
    
    # shift to positive
    corF <- apply(corF, 2, function(x){ x + abs(min(x)) })
    
    # fill the function matrix
    Funcs[1:nrow(Funcs), 1:ncol(Funcs)] <- 
      as.vector(corF)
    
    # convert the Funcs data to a dataframe
    Funcs <- as.data.frame(Funcs)
    
    # add a species column
    Funcs$species <- row.names(Funcs)
    row.names(Funcs) <- NULL
    
    # add these function values to the lotka_volterra data
    df <- 
      left_join(x,
                Funcs,
                by = "species") %>%
      mutate( across(.cols = all_of(func_names), ~(.*abundance) ) ) %>%
      group_by(replicate, species_pool) %>%
      summarise( across(.cols = c("abundance", all_of(func_names) ) , ~sum(.) ), .groups = "drop" )
    
    # z-score standardise the different functions
    # translate them to make sure they are positive
    df <- 
      df %>%
      mutate(across(.cols = all_of(func_names), ~as.numeric(scale(.x, center = TRUE, scale = TRUE)) )) %>%
      mutate(across(.cols = all_of(func_names), ~(.x + abs(min(.x)))  ))
    
    ss_funcs[[i]] <- df
    
  }
  
  bind_rows(ss_funcs, .id = "cor_mat")
  
})

# bind this into a dataframe
lv_mf <- bind_rows(lv_mf, .id = "model")

# write this into a .csv file
write_csv(x = lv_mf,
          path = here("data/lv_mf_analysis_data.csv"))






