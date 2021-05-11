
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


# set up different functional scenarios for different species

# set the number of functions
n_funcs <- 9

# set the function names
func_names <- paste0("F_", 1:n_funcs)

# scenario 1:
# "all species positively affect all functions as a function of their abundance"

# scenario 2:
# "species have a mix of positive and negative effects on different functions as a function of their abundance"
# "however, positive effects outweigh negative effects"

# scenario 3:
# "species have a mix of positive and negatives on different functions as a function of their abundance

# defines lower bound of the uniform distribution for scenario 1, 2 and 3 respectively
mixf <- c(0, -0.25, -0.5)

sp_funcs <- 
  lapply(mixf, function(x) {
    
    z <- 
      data.frame(t(replicate(n = length(sp_names), expr = runif(n = n_funcs, min = x, max = 1))))
    
    names(z) <- func_names
    
    z$species <- sp_names
    
    z
    
  })

# for each of these species-specific function values, multiply it by the species' abundances
# sum all species values in each community
# standardise by the max
com_mf <- 
  lapply(sp_funcs, function(x) {
    
    y <- left_join(lv_dat, x, by = "species")
    
    y <- 
      y %>%
      mutate(across(.cols = all_of(func_names), ~(.*abundance) ))
    
    y <- 
      y %>%
      group_by(run, species_pool, replicate) %>%
      summarise( across(.cols = c("abundance", all_of(func_names) ) , 
                        ~if_else(sum(.) >= 0, sum(.), 0)),
                 .groups = "drop" )
    
    y <- 
      y %>%
      group_by(run) %>%
      mutate( across(.cols = all_of(func_names) , ~./max(.)) ) %>%
      ungroup()
    
    y
    
  })


# explore the com_mf data
com_mf

# plot the relationship between species pool diversity and function
com_mf[[1]] %>%
  pivot_longer(cols = all_of(func_names),
               names_to = "eco_function",
               values_to = "function_value") %>%
  mutate(run = as.character(run)) %>%
  ggplot(data = .,
       mapping = aes(x = species_pool, y = function_value, colour = run)) +
  geom_jitter(alpha = 0.1, width = 0.2) +
  geom_smooth(se = FALSE, method = "lm") +
  facet_wrap(~eco_function, scales = "free") +
  theme_meta()

com_mf[[2]] %>%
  filter(run == 3) %>%
  select(all_of(func_names)) %>%
  cor() %>%
  corrplot()

# prepare the data for export
names(com_mf) <- paste0("s_", 1:length(mixf))

com_mf <- bind_rows(com_mf, .id = "scenario")

# write this into a .csv file
write_csv(x = com_mf,
          path = here("data/lv_mf_analysis_data.csv"))






