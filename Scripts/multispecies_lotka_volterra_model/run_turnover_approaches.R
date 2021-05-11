
# Project: Review of multifunctionality in ecology, conservation and ecosystem service science

# Title: Prepare cluster of model replicates (i.e. 1000 simulations of one set of parameters)

# load relevant libraries
library(here)

# link to scripts with the relevant functions
source(here("Scripts/multispecies_lotka_volterra_model/process_model_data.R"))
source(here("Scripts/multispecies_lotka_volterra_model/per_capita_function_matrix_functions.R"))
source(here("Scripts/multispecies_lotka_volterra_model/turnover_approach_functions.R"))

# Laura will provide function matrices so that the identity effects in the DI-models are the same
# for now, I will use the func_matrix_generator function:

# choose the number of functions
n.f <- 5

# specialist
set.seed(123)
fm <- func_matrix_generator(species_list = unique(x$species), 
                            func.n = n.f, func.spec = "specialist", prob.neg = 0.05)


# simulated data cluster (e.g. neutral model with 1000 runs with same parameters)
x

# process the simulated cluster using: process_sim_data
data_proc <- process_sim_data(model_data = x, func.mat = fm, time_final = TRUE, species_abun = "pa")


# function to compare AIC and SES-based turnover approaches: turnover_tester()

# arguments:
# model_dat: data.frame from functions: drift_model.R or stachova_leps_2010
# that has been processed using the process_sim_data() function
# function_matrix: data.frame with identity effect coefficients for all species in model_dat

turnover_tester <- 
  function(model_dat, function_matrix, ses_reps = 100) {
  
  # check that the correct packages are installed
  if(! "dplyr" %in% installed.packages()[,1]) stop(
    "this function requires dplyr to be installed"
  )
    
  # load the dplyr library
  library(dplyr)
    
  # convert this data cluster into a list
  model.list <- split(model_dat, model_dat$model_run)
  
  # get a vector of function names from the function matrix
  f.names <- names(function_matrix)[grepl(pattern = "F_", names(function_matrix))]
  
  turnover_list <- 
    lapply(model.list, function(data.mf) {
      
      # get a vector of species names that are present for the chosen simulation
      sp.present <- sapply(data.mf[, grepl("sp_", names(data.mf)) ], function(x) sum(ifelse(x > 0, 1, 0)))
      sp.present <- names(sp.present[sp.present > 0])
      
      # subset the species that are present in the simulation
      data.in <- 
        data.mf %>%
        select(model_run, patch, time, local_species_pool, composition,
               all_of(sp.present), all_of(f.names))
      
      # subset the present species
      func.in <- fm[fm$species %in% sp.present, ]
      
      
      ### AIC-based turnover approach
      
      # implement the aic-based approach to get species effects on each function
      aic.x <- AIC_sp(data = data.in, function_names = f.names, species_names = sp.present)
      
      # calculate the proportion of incorrect directions for each function
      aic.dr <- mapply(compare.directions, func.in[, f.names], aic.x[, f.names])
      
      # calculate the spearman correlation
      aic.spear <- mapply(function(x, y){cor(x, y, method = "spearman")}, func.in[, f.names], aic.x[, f.names])
      
      
      ### SES-based turnover approach
      
      # implement the ses-based approach to get species effects on each function
      ses.x <- SES_score(data = data.in, function_names = f.names, species_names = sp.present, n_ran = ses_reps )
      
      # calculate the proportion of incorrect directions for each function
      ses.dr <- mapply(compare.directions, func.in[, f.names], ses.x[, f.names])
      
      # calculate the spearman correlation
      ses.spear <- mapply(function(x, y){cor(x, y, method = "spearman")}, func.in[, f.names], ses.x[, f.names])
      
      # summary metric
      summary_stats <- 
        bind_rows(aic.dr, aic.spear, ses.dr, ses.spear) %>%
        mutate(output_metric = rep(c("prop_incorrect", "spearman_r"), 2),
               method = rep(c("AIC", "SES"), each = 2)) %>%
        select(method, output_metric, all_of(f.names))
      
      return(summary_stats)
      
    })
  
  # bind the output list into a dataframe
  turnover_summary <- bind_rows(turnover_list, .id = "model_run")
  
  return(turnover_summary)
  
}

# test the turnover_tester() function
turnover_tester(model_dat = mf.turnover, function_matrix = fm)

### END
