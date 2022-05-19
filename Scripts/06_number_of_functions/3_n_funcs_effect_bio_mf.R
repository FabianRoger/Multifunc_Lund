
# Project: Review of multifunctionality in ecology, conservation and ecosystem service science

# Title: Does the number of functions matter?

# load relevant libraries
library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)
library(here)

rm(list = ls())

# load functions from important scripts
source(here("Scripts/01_general_functions/MF_functions_collated.R"))
source(here("Scripts/01_general_functions/process_model_data.R"))


# using functions defined to get combinations of functions:

# set up a function to calculate the realised diversity-multifunctionality slope
# for each combination of a set of ecosystem functions

# adf.data: data.frame where rows are plots and which contains ecosystem functions of interest
# mf.func.names: vector of function names to consider
# standardise_funcs = TRUE i.e. standardise functions (z-score) before calculating multifunctionality
# covariate_name = name of the variable to regress against the multifunctional BEF slope

get_BEF_mf_est <- function(adf.data, 
                           mf.func.names, 
                           standardise_funcs = TRUE,
                           covariate_name = "realised_diversity") {
  
  if (length(mf.func.names) < 2) {
    stop("must have more than two functions to adequately examine the effect of changing the number of functions")
  }
  
  # write the adf.data into a new data.frame
  dat.in <- adf.data
  
  # use the function to get a list of all function combinations
  list.func.names <- get.function.combinations(function.names = mf.func.names)
  
  # set an output list
  list.out <- vector("list", length(list.func.names))
  
  # loop over each set of function names
  for (i in 1:length(list.func.names)) {
    
    # get the vector of function names
    sample.func.names <- list.func.names[[i]]
    
    # calculate the multifunctionality metrics
    data.mf <- 
      dat.in %>%
      mutate(sum_MF = MF_sum(adf = dat.in, vars = sample.func.names),
             ave_MF = MF_av(adf = dat.in, vars = sample.func.names, stand_method = "z_score_abs"),
             sd_MF = MF_sd(adf = dat.in, vars = sample.func.names, stand_method = "z_score_abs"),
             Pasari_MF = MF_pasari(adf = dat.in, vars = sample.func.names),
             ENF_MF = hill_multifunc(adf = dat.in, vars = sample.func.names, scale = 1, HILL = TRUE),
             thresh_30_MF = single_threshold_mf(adf = dat.in, vars = sample.func.names, thresh = 0.3),
             thresh_70_MF  = single_threshold_mf(adf = dat.in, vars = sample.func.names, thresh = 0.7) 
             )
    
    # calculate summary statistics for the data.mf data.frame
    data.mf.summary <- 
      data.mf %>%
      pivot_longer(cols = ends_with("MF"),
                   names_to = "multifunctionality_metric",
                   values_to = "MF") %>%
      dplyr::select(multifunctionality_metric, MF) %>%
      group_by(multifunctionality_metric) %>%
      summarise(mean_MF = mean(MF, na.rm = TRUE),
                sum_MF = sum(MF, na.rm = TRUE),
                sd_MF = sd(MF, na.rm = TRUE),
                cv_MF = (sd_MF/mean_MF),
                range_MF = diff(range(MF)),
                min_MF = min(MF, na.rm = TRUE),
                max_MF = max(MF, na.rm = TRUE), .groups = "drop")
    
    # subset the multifunctionality metric names
    mf.names <- names(data.mf)[grepl(pattern = "MF", x = names(data.mf))]
    
    # calculate summary statistics for the among function variation
    data.f.summary <- 
      dat.in %>%
      mutate(row_id = 1:nrow(dat.in)) %>%
      dplyr::select(row_id, local_species_pool, all_of(sample.func.names)) %>%
      pivot_longer(cols = all_of(sample.func.names),
                   names_to = "function_id",
                   values_to = "function_value") %>%
      group_by(local_species_pool, row_id) %>%
      summarise(F_mean = mean(function_value, na.rm = TRUE),
                F_sd = sd(function_value, na.rm = TRUE),
                F_cv = (sd(function_value, na.rm = TRUE)/mean(function_value, na.rm = TRUE)),
                F_range = diff(range(function_value)), .groups = "drop") %>%
      summarise(sd_funcs = mean(F_sd),
                cv_funcs = mean(F_cv),
                range_funcs = mean(F_range),
                mean_cor = cor(local_species_pool, F_mean, method = "spearman"),
                sd_cor = cor(local_species_pool, F_sd, method = "spearman"),
                cv_cor = cor(local_species_pool, F_cv, method = "spearman"),
                range_cor = cor(local_species_pool, F_range, method = "spearman"))
    
    # for each multifunctionality metric, calculate the BEF-slope
    bef_mf_slope <- 
      sapply(mf.names, function(x){
        lm.x <- lm(reformulate(covariate_name, as.name(x) ), data = data.mf)
        lm.x$coefficients[2]
      })
    
    names(bef_mf_slope) <- NULL # remove names from the vector
    
    # bind this into a data.frame
    df.out <- 
      data.frame(func.comb.id = i,
                 n.func.id = paste(sample.func.names, collapse = ""),
                 number_of_functions = length(sample.func.names),
                 sd_funcs = data.f.summary$sd_funcs,
                 cv_funcs = data.f.summary$cv_funcs,
                 range_funcs = data.f.summary$range_funcs,
                 mean_cor = data.f.summary$mean_cor,
                 sd_cor = data.f.summary$sd_cor,
                 cv_cor = data.f.summary$cv_cor,
                 range_cor = data.f.summary$range_cor,
                 multifunctionality_metric = mf.names,
                 diversity_mf_est = bef_mf_slope)
    
    # join the summary statistics to the BEF slope data
    df.out <- full_join(df.out, data.mf.summary, by = "multifunctionality_metric")
    
    list.out[[i]] <- df.out
  }
  
  bind_rows(list.out)
  
}

# load the processed model data
mod.out <- read_csv(here("data/drift_model_n_functions_processed.csv"))
head(mod.out)
unique(mod.out$function_matrix)
length(unique(mod.out$mod_id))

# do this for a subset of the data i.e. n = 50
mod.out <- 
  mod.out %>%
  filter(mod_id %in% sample(unique(mod.out$mod_id), 200))

# output the function names
f.names <- names(mod.out)[grepl(pattern = "F_", names(mod.out))]

# get the BEF-multifunctional slope for each simulated dataset across all metrics
# note this code can take up a long time to run (e.g. 20-40 mins depending on the speed of your computer)

sim.n.func <- 
  lapply(split(mod.out, mod.out$mod_id), function(data) {
  
  get_BEF_mf_est(adf.data = data, 
                 mf.func.names = f.names, 
                 standardise_funcs = TRUE,
                 covariate_name = "local_species_pool")
  
})

sim.n.df <- bind_rows(sim.n.func, .id = "mod_id")
head(sim.n.df)
names(sim.n.df)
View(sim.n.df)

# join the identifying variables
sim.n.out <- 
  left_join(sim.n.df,
            mod.out %>%
              dplyr::select(mod_id, drift_parameter, function_matrix, model_run) %>%
              distinct(), by = "mod_id")
head(sim.n.out)

# write this to a csv so we don't have to run the model again
write_csv(x = sim.n.out, file = here("data/sim_n_functions.csv"))

### END
