
# Project: Examining the relationship between biodiversity and ecosystem functioning in experimental and observational data

# Title: Preparation of the stachova and leps model output data for the multifunctionality work

# load relevant libraries generally
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(viridis)
library(here)

# where to access functions from
source(here("Scripts/multispecies_lotka_volterra_model/multispecies_lotka_volterra_model.R"))

# run the model to generate the output data

# set the number of runs to do
n_exp <- 5

# set up the means for the model runs
set.seed(4897245)
a_mean_sim <- runif(n = n_exp, min = 0.25, max = 1)

sl_mod_out <- vector("list", length = n_exp)
for (i in (1:n_exp) ) {
  
  sl_mod_out[[i]] <- 
    s_l_2010_mod(reg_pool = 20,
                 t_steps = 4000, 
                 n0 = 3,
                 a_mean = a_mean_sim[i], 
                 a_sd = 0.2, a_min = 0.2, a_max = 1.2, a_spp = 1,
                 k_min = 3, k_max = 150,
                 r_min = 0.01, r_max = 0.5, 
                 lsp = c(1, 2, 4, 6, 8, 10, 12),
                 reps = 20)
  
}

sl_mod_an <- 
  bind_rows(sl_mod_out, .id = "run")

# filter out the final time point in the model
sl_mod_an <- 
  sl_mod_an %>%
  filter(time == last(time))

# output this into dataframe as a .csv file
write_csv(x = sl_mod_an,
          path = here("data/stachova_leps_model_mf.csv"))

