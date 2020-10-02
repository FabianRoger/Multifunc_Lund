
# Project: Mutifunctionality workshop (Lund 2019)

# Title: Adding multifunctionality into Lotka-Volterra models

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

# function to get the slope for different numbers of functions in a dataset
# data is a dataset with multiple functions as columns
# funcs is a vector of function names
est_n_func <- function(data, funcs) {
  
  s_est <- function(a, b, d) {
    
    # set up combinations for each number of functions
    comb_f <- gtools::combinations(n = length(b), r = d)
    
    # set up an output vector
    y <- vector("list", length = nrow(comb_f))
    
    for (i in 1:nrow(comb_f)) {
      
      v <- b[comb_f[i, ]]
      
      df <- a[, v]
      
      lm_df <- 
        a %>%
        mutate(species_pool = scale(species_pool, scale = TRUE, center = TRUE),
               ave_mf = MF_av(adf = df, vars = v),
               sum_mf = MF_sum(adf = df, vars = v),
               pasari_mf = MF_pasari(adf = df, vars = v),
               enf_mf = hill_multifunc(adf = df, vars = v, scale = 1, HILL = TRUE),
               thresh_30_mf = single_threshold_mf(adf = df, vars = v, thresh = 0.3),
               thresh_70_mf = single_threshold_mf(adf = df, vars = v, thresh = 0.7))
      
      mf_names <- names(select(lm_df, ends_with("mf")))
      
      s_cof <- vector(length = length(mf_names))
      
      for(s in 1:length(mf_names)) {
        
        lm1 <- lm(as.formula(paste0(mf_names[s], "~ species_pool")),
                  data = lm_df)
        
        s_cof[s] <- lm1$coefficients[2]
        
      }
      
      names(s_cof) <- c(mf_names)
      
      y[[i]] <- s_cof
      
    } 
    
    y <- bind_rows(y, .id = "run")
    
    y$number_functions <- d
    
    y
    
  }
  
  # loop over this function for all possible numbers of functions
  slopes_n_func <- vector("list", length = length(funcs))
  
  for(j in 2:length(funcs)) {
    
    slopes_n_func[[j]] <- 
      s_est(a = data, b = funcs, d = j)
    
  }
  
  bind_rows(slopes_n_func, .id = "n_func")
  
}


# load the Lotka-Volterra simulated data
lv_dat <- 
  read_csv(file = here("data/stachova_leps_model_mf.csv"))

# check dimensions of the data
nrow(lv_dat)

# make a vector of species names
sp_names <- unique(lv_dat$species)

# choose the number of functions to simulate
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
  
  # set correlation cluster values
  lc <- 
    list(c1 = c(0, 0, 0),
         c2 = c(-0.4, -0.3, -0.5),
         c3 = c(0.1, 0.4, 0.4))
  
  
  # make three 'cluster' of correlated functions
  Sigma[1:3,1:3] <- lc[[2]][1]
  Sigma[8:9,8:9] <- lc[[2]][2]
  Sigma[6:7,6:7] <- lc[[2]][3]
  
  diag(Sigma) <- 1
  
  is.positive.definite(x = Sigma)
  
  # draw correlated functions (with mean 0)
  corF <- mvrnorm(n = specnum, mu = rep(0, func_n), Sigma = Sigma)
  
  # shift to positive
  corF <- apply(corF, 2, function(x){ x + abs(min(x)) })
  
  # fill the function matrix
  Funcs[1:nrow(Funcs), 1:ncol(Funcs)] <- 
    as.vector(corF)
  
  Funcs %>%
    cor() %>%
    corrplot(method = "ellipse")
  
  # convert the Funcs data to a dataframe
  Funcs <- as.data.frame(Funcs)
  
  # add a species column
  Funcs$species <- row.names(Funcs)
  row.names(Funcs) <- NULL
  
  
  # extract the abundance of each species in the last time point
  # multiply each species abundance by the function value
  # calculate function scores for each species
  
  bio_funcs <- 
    left_join(filter(lv_dat, run == 1),
              Funcs,
              by = "species") %>%
    mutate( across(.cols = all_of(func_names), ~(.*abundance) ) ) %>%
    group_by(replicate, species_pool) %>%
    summarise( across(.cols = c("abundance", all_of(func_names) ) , ~sum(.) ), .groups = "drop" )
  
  # z-score standardise the different functions
  # translate them to make sure they are positive
  bio_funcs <- 
    bio_funcs %>%
    mutate(across(.cols = all_of(func_names), ~as.numeric(scale(.x, center = TRUE, scale = TRUE)) )) %>%
    mutate(across(.cols = all_of(func_names), ~(.x + abs(min(.x)))  ))
  
  bio_funcs %>%
    select(abundance, contains("F")) %>%
    cor() %>%
    corrplot(method = "ellipse")
  
  # plot the diversity-function relationship
  bio_funcs %>%
    pivot_longer(cols = starts_with("F"),
                 names_to = "eco_function",
                 values_to = "value") %>%
    ggplot(data = .,
           mapping = aes(x = species_pool, y = value)) +
    geom_point() +
    geom_smooth(method = "lm", se = FALSE) +
    facet_wrap(~eco_function) +
    theme_meta()
  
  # run the est_n_function to get the slope between multifunctionality and richness
  df_est <- est_n_func(data = bio_funcs, funcs = func_names)
  
  # plot the standardised effect size between species pool diversity and function
  df_est <- 
    df_est %>%
    pivot_longer(cols = ends_with("mf"),
                 names_to = "mf_metric",
                 values_to = "biodiv_func_est")
  
ggplot(data = df_est, 
       mapping = aes(x = number_functions, biodiv_func_est)) +
  geom_jitter(width = 0.5, alpha = 0.1) +
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(~mf_metric, scales = "free") +
  theme_meta()

ggplot(data = df_est %>%
         filter(grepl("thresh", mf_metric)), 
       mapping = aes(x = number_functions, biodiv_func_est)) +
  geom_jitter(width = 0.5, alpha = 0.1) +
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(~mf_metric) +
  theme_meta()






















