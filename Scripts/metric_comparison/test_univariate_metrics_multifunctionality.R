
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
  
}

# simulate data for five ecosystem functions
sim.dat <- univariate_explorer() 



