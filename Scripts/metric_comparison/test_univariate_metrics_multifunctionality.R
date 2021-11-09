
# Project: Review of multifunctionality in ecology, conservation and ecosystem service science

# Title: Code to compare different metrics of multifunctionality given different function distributions

# load relevant libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(here)

rm(list = ls())

# tell R where to call scripts from
source(here("Scripts/function_plotting_theme.R"))
source(here("Scripts/MF_functions_collated.R"))

# write a function to generate these function distributions

univariate_explorer <- function(funcnum = 5, grain = 0.05, error = 0.01) {
  
  # generate n functions between 0 and 1
  x <- replicate(funcnum, seq(0, 1, grain), simplify = FALSE)
  
  # put these functions into a matrix
  y <- do.call("expand.grid", x)
  
  # add a small amount of error and set to zero if zero
  nm <- nrow(y)*ncol(y)
  y <- y + rnorm(n = nm, 0, error)
  y[y < 0] <- 0
  
  names(y) <- paste("F_", 1:funcnum, sep = "")
  
  # convert y to a tibble
  df <- tibble(y)
  
  # only get unique rows
  df <- distinct(df)
  
  # output the df
  df
  
}

# simulate data for five ecosystem functions
sim.dat <- univariate_explorer() 
sim.dat

# get a vector of names
func.names <- names(sim.dat)

sim.dat <- 
  sim.dat %>%
  mutate(`scal. MF` = MF_jing(adf = sim.dat, vars = func.names),
         `sum MF` = MF_sum(adf = sim.dat, vars = func.names),
         `ave. MF` = MF_av(adf = sim.dat, vars = func.names),
         `MESLI MF` = MF_mesli(adf = sim.dat, vars = func.names),
         `Pasari MF` = MF_pasari(adf = sim.dat, vars = func.names),
         `SAM MF` = MF_dooley(adf = sim.dat, vars = func.names),
         `ENF MF` = hill_multifunc(adf = sim.dat, vars = func.names, scale = 1, HILL = TRUE),
         `Simp. MF` = MF_simpsons_div(adf = sim.dat, vars = func.names),
         `Manning.30 MF` = manning_multifunc(adf = sim.dat, vars = func.names, thresh = 0.3),
         `Manning.50 MF` = manning_multifunc(adf = sim.dat, vars = func.names, thresh = 0.5),
         `Manning.70 MF` = manning_multifunc(adf = sim.dat, vars = func.names, thresh = 0.7),
         `thresh.30 MF` = single_threshold_mf(adf = sim.dat, vars = func.names, thresh = 0.3),
         `thresh.50 MF` = single_threshold_mf(adf = sim.dat, vars = func.names, thresh = 0.5),
         `thresh.70 MF` = single_threshold_mf(adf = sim.dat, vars = func.names, thresh = 0.7),
         `Slade.10.90 MF` = MF_slade(adf = sim.dat, vars = func.names, A_quant = 0.10, B_quant = 0.90),
         `Slade.40.60 MF` = MF_slade(adf = sim.dat, vars = func.names, A_quant = 0.40, B_quant = 0.60),
         `PCA MF` = pca_multifunc(adf = sim.dat, vars = func.names, standardise = FALSE) 
         )

# save this as a .rds file
saveRDS(sim.dat, here("Scripts/metric_comparison/metric_sims.rds") )







