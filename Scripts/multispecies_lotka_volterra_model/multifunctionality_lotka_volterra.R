
# Project: Mutifunctionality workshop (Lund 2019)

# Title: Adding multifunctionality into Lotka-Volterra models

# load libraries
library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)
library(here)
library(corrplot)

# load the Lotka-Volterra simulated data
lv_dat <- 
  read_csv(file = here("Scripts/multispecies_lotka_volterra_model/lv_mf_sim_dat.csv"))

# check dimensions of the data
nrow(lv_dat)

# make a vector of species names
sp_names <- unique(lv_dat$species)

# choose the number of functions to simulate
func_n <- 5

# set the function names
func_names <- paste0("F_", 1:func_n)

# for each species, assign a function value for x functions

# create an empty data matrix
Funcs <- 
  matrix(
  nrow = length(sp_names),
  ncol = func_n,
  dimnames = list(sp_names, func_names))

# fill this empty data matrix with values from the uniform distribution
Funcs[1:nrow(Funcs), 1:ncol(Funcs)] <- 
  runif( n = prod( dim( Funcs ) ), 0, 1 )

# convert the Funcs data to a dataframe
Funcs <- as.data.frame(Funcs)

# add a species column
Funcs$species <- row.names(Funcs)
row.names(Funcs) <- NULL


# extract the abundance of each species in the last time point
# multiply each species abundance by the function value
# calculate function scores for each species

bio_funcs <- 
  left_join(filter(lv_dat, time == last(time)),
            Funcs,
            by = "species") %>%
  mutate( across(.cols = func_names, ~(.*abundance) ) ) %>%
  group_by(replicate, species_pool) %>%
  summarise( across(.cols = c("abundance", all_of(func_names) ) , ~sum(.) ), .groups = "drop" )

bio_funcs %>%
  select(abundance, contains("F")) %>%
  cor() %>%
  corrplot(method = "ellipse")

















