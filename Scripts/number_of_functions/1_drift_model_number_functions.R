
# Project: Review of multifunctionality in ecology, conservation and ecosystem service science

# Title: Run the neutral model to test the number of functions question

# load relevant libraries
library(here)

# clear the current memory
rm(list = ls())

# link to scripts with the relevant functions
source(here("Scripts//ecological_drift_model.R"))
source(here("Scripts/MF_functions_collated.R"))

# how many model reps for each parameter combination?
# with 1000 replicates, it can take up to 30 minutes to complete
n_reps <- 1000

# run the ecological drift model

# drift model parameters
p_change <- 0.025

drift.mod.list <- vector("list", length = n_reps)
for (i in 1:n_reps) {
  
  library(dplyr)
  
  drift.mod.list[[i]] <- 
  
    drift_model(lsp = c(2, 4, 6, 9),
                mono = "all",
                reps = 5,
                technical_reps = 2,
                rsp = 12,
                t_steps = 500,
                n0 = 500,
                prop_change = p_change,
                n_repeats = 1) %>%
    filter(time == max(time))
}

# bind the list into a data.frame
drift.mod.list <- bind_rows(drift.mod.list, .id = "model_run")
drift.mod.list$drift_parameter <- p_change

# reorganise the columns
drift.mod.list <- 
  drift.mod.list %>%
  select(drift_parameter, model_run, time, patch, local_species_pool, composition,
         species, abundance)

# write this model data.frame into a .csv file so the model does not have to be re-run
library(readr)

# make a folder to export the cleaned data
if(! dir.exists(here("data"))){
  dir.create(here("data"))
}

write_csv(x = drift.mod.list, file = here("data/drift_model_n_functions.csv"))

### END

