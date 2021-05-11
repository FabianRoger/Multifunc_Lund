
# Project: Review of multifunctionality in ecology, conservation and ecosystem service science

# Title: Function to process the raw model data for the turnover approach

# arguments:
# model_data: cluster of data outputted from drift_model() or s_l_function()
# func.mat: function matrix outputed from per_capita_function script
# time_final: if TRUE then only the final time point gets outputted otherwise all time-points are retained
# species_abun:
# "pa": converts species abundances to presence-absences (0, 1)
# "proportion": converts species abundances to proportions within each patch (0 to 1)
# "raw": keeps the raw abundance data i.e. the counts

# get the simulated data cluster (i.e. n reps of a model with the same parameters)
process_sim_data <- function(model_data, func.mat, time_final = TRUE, species_abun = "pa") {
  
  # check that the correct packages are installed
  if(! "dplyr" %in% installed.packages()[,1]) stop(
    "this function requires dplyr to be installed"
  )
  
  # check that the correct packages are installed
  if(! "tidyr" %in% installed.packages()[,1]) stop(
    "this function requires tidyr to be installed"
  )
  
  # load the dplyr library
  library(dplyr)
  library(tidyr)
  
  # write the model data set into an object called df 
  df <- model_data
  
  # write the function matrix into an object called fm
  fm <- func.mat
  
  # extract a vector of function names from fm
  f.names <- names(fm[,-1])
  
  # extract the final time-point: TRUE or FALSE
  if(time_final == TRUE) {
    
    df <- 
      df %>%
      filter(time == max(time))
    
  }
  
  
  # (1) join the abundance data to the function values
  # (2) multiply function coefficients by species abundances
  # (3) add normally distributed error the per species function values
  df_mf <- 
    full_join(df, fm, by = "species") %>%
    dplyr::mutate( dplyr::across(.cols = all_of(f.names) , ~(.*abundance) ) ) %>%
    dplyr::mutate( dplyr::across(.cols = all_of(f.names) , ~(. + rnorm(n = length(.), mean = 0, sd = 0.1)) ) )
  
  # summarise the data to get per patc values of abundance and 
  df_mf_s <- 
    df_mf %>%
    group_by(model_run, patch, time, local_species_pool, composition) %>%
    summarise_at(vars( all_of(c("abundance", f.names)) ), sum) %>%
    ungroup()
  
  # convert these data into the wide format
  df_mf_w <- 
    df_mf %>%
    select(-all_of(f.names))
  
  # transform the abundance data into either presence-absence of proportions
  if (species_abun == "pa") {
    
    df_mf_w <- 
      df_mf_w %>%
      mutate(abundance = if_else(abundance > 0, 1, 0))
    
  } else if (species_abun == "proportion") {
    
    df_mf_w <- 
      df_mf_w %>%
      group_by(model_run, patch, time, local_species_pool, composition) %>%
      mutate(abundance = abundance/sum(abundance)) %>%
      ungroup()
    
  } else if (species_abun == "raw") {
    
  } else { 
    
    stop("error! choose appropriate species_abun argument") 
    
  }
  
  # convert these data into the wide format
  df_mf_w <- 
    df_mf_w %>%
    pivot_wider(id_cols = everything(),
                names_from = "species", values_from = "abundance")
  
  # join these datasets together
  mf_data <- full_join(df_mf_w, df_mf_s, by = c("model_run", "patch", "time", "local_species_pool", "composition"))
  
  return(mf_data)
  
}

### END
