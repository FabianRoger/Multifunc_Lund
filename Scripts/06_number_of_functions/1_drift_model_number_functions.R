
# Project: Review of multifunctionality in ecology, conservation and ecosystem service science

# Title: Run the neutral model to test the number of functions question

# load relevant libraries
library(here)
library(readr)

# clear the current memory
rm(list = ls())

# link to scripts with the relevant functions
source(here("Scripts/01_general_functions/ecological_drift_model.R"))
source(here("Scripts/01_general_functions/MF_functions_collated.R"))

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
                n0 = 1000,
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

# make a folder to export the cleaned data
if(! dir.exists(here("data"))){
  dir.create(here("data"))
}

# write_csv(x = drift.mod.list, file = here("data/drift_model_n_functions.csv"))


# Test the drift model i.e. is there a relationship between abundance and species richness
drift.mod.list <- read_csv(here("data/drift_model_n_functions.csv"))

# calculate the slope between number of functions and multifuncitonal BEF slope 
# and other summary statistics for each multifunctionality metric
library(purrr)
library(broom)
library(tidyr)
library(ggplot2)

# define function to efficiently output the slope
lm.cleaner <- function(data, response, explanatory, output_prefix = "x") {
  
  x <- 
    lm(reformulate(explanatory, response), data = data) %>% 
    tidy %>% 
    filter(term == explanatory) %>% 
    select(!!paste(output_prefix, response, sep = "") := estimate )
  
  return(x)
  
}

df_test <- 
  drift.mod.list %>%
  group_by(model_run, local_species_pool, patch) %>%
  summarise(abundance = sum(abundance)) %>%
  group_by(model_run) %>% 
  mutate(across(.cols = c("abundance"), ~standardise_functions(x = ., method = "z_score") )) %>%
  nest() %>% 
  summarise(total_abun_slope = map(data, ~lm.cleaner(data = .x, explanatory = "local_species_pool", response = "abundance")) ) %>%
  unnest(ends_with("slope"))
head(df_test)

# check the range of slopes
range(df_test$xabundance)

# plot a histogram of slopes
ggplot(data = df_test,
       mapping = aes(x = xabundance)) +
  geom_histogram(colour = "transparent", alpha = 0.5) +
  geom_vline(xintercept = mean(df_test$xabundance), colour = "red") +
  geom_vline(xintercept = 0, colour = "black", linetype = "dashed") +
  theme_classic()

### END
