
# Project: Review of multifunctionality in ecology, conservation and ecosystem service science

# Title: Prepare cluster of model replicates (i.e. 1000 simulations of one set of parameters)

# next steps:
# add packages to 'process_model_data'
# test the functions and see if they can consistently output the relevant metrics

# get Laura to send me a functional matrix and the data in this format form the DI models

# load relevant libraries
library(dplyr)
library(tidyr)

# link to scripts with the relevant functions
# source() etc.
# compare.directions
# AIC_sp
# SES_score
# func_matrix_generator
# process_sim_data

# simulate specialist and generalist function matrices

# choose the number of functions
n.f <- 5

# specialist
set.seed(123)
fm <- func_matrix_generator(species_list = unique(x$species), 
                            func.n = n.f, func.spec = "specialist", prob.neg = 0.05)

# generalist
set.seed(123)
fm <- func_matrix_generator(species_list = unique(x$species), 
                            func.n = n.f, func.spec = "generalist", prob.neg = 0.05)


# get the simulated data cluster (i.e. n reps of a model with the same parameters)
x

# process the simulated cluster using: process_sim_data
process_sim_data(model_data = x, func.mat = fm, time_final = TRUE, species_abun = "pa")

# convert this data cluster into a list
mf.turnover.list <- split(mf.turnover, mf.turnover$model_run)

l.out <- 
  lapply(mf.turnover.list, function(data.mf) {
  
  # get a vector of species names that are present for the chosen simulation
  sp.present <- sapply(data.mf[, grepl("sp_", names(data.mf)) ], function(x) sum(ifelse(x > 0, 1, 0)))
  sp.present <- names(sp.present[sp.present > 0])
  
  # subset the species that are present in the simulation
  data.in <- 
    data.mf %>%
    select(model_run, patch, time, local_species_pool, composition,
           all_of(sp.present), all_of(f.names))
  
  # subset the present species
  func.in <- fm[fm$species %in% sp.present,]
  
  # implement the aic-based approach to get species effects on each function
  aic.x <- AIC_sp(data = data.in, function_names = f.names, species_names = sp.present)
  
  aic.test <- mapply(compare.directions, func.in[, f.names], aic.x[, f.names])
  
  # implement the ses-based approach to get species effects on each function
  ses.x <- SES_score(data = data.in, function_names = f.names, species_names = sp.present, n_ran = 100)
  
  ses.test <- mapply(compare.directions, func.in[, f.names], ses.x[, f.names])
  
  test.bind <- 
    bind_rows(aic.test, ses.test) %>%
    mutate(method = c("AIC", "SES")) %>%
    select(method, all_of(f.names))
  
  return(test.bind)
  
})

bind_rows(l.out, .id = "model_run")


