
# Project: Review of multifunctionality in ecology, conservation and ecosystem service science

# Title: Implement the turnover approach using the multifunc package and an analogue to Meyer et al. (2017)'s null model

# set-up general functions

# function to mix-up relationship between function and species
row.randomiser <- function(func.names, species.names, adf.data) {
  
  # subset a matrix of functions
  func.mat <- adf.data[, func.names]
  func.mat.nrow <- nrow(func.mat) # calculat the number of rows in that matrix
  
  # use sample to get random row ids from the function matrix
  random.row.ids <- sample(x = 1:func.mat.nrow, size = func.mat.nrow , replace = FALSE)
  
  # randomise the function matrix
  func.mat.random <- func.mat[random.row.ids, ]
  
  # get a matrix of species data
  spec.mat <- adf.data[, species.names]
  
  # bind this randomised data together
  adf.data.random <- cbind(spec.mat, func.mat.random)
  
  return(as.data.frame(adf.data.random))
  
}

# function to test reliability

# Q1: for detected species effects, how many species effects are in the incorrect direction?

# true.function.vals: real function values (i.e. linear coefficient)
# inferred.effects: inferred species effect (-1, 0 or 1)

# note that this function ignores non-significant species effects (i.e. effect = 0)

# outputs proportion of IN-correctly inferred directions
compare.directions <- function(true.function.vals, inferred.effects){
  n.incorrect.negative <- sum(true.function.vals[inferred.effects < 0] > 0)
  n.incorrect.positive <- sum(true.function.vals[inferred.effects > 0] < 0)
  n.incorrect.total <- (n.incorrect.negative + n.incorrect.positive)
  
  return( n.incorrect.total/sum(inferred.effects != 0) )
}


# functions to implement aic-based species effects and Gotelli et al. (2011) species effects

# 1. calculate species effects using the AIC-based turnover approach

# data: data.frame with functions and either species abundances or presence-absences as columns and samples as rows
# function_names: vector of function names
# species_names: vector of species names

AIC_sp <- function(data, function_names, species_names) {
  
  if(! "multifunc" %in% installed.packages()[,1]) stop(
    "this function requires the multifunc package to be installed"
  )
  
  # load the multifunc package
  library(multifunc)
  
  sp.effect.aic <- vector("list", length = length(function_names))
  for (i in 1:length(function_names)){
    
    # get species effect on each function using abundances
    # this outputs a vector of species effects (-1, 0 or 1) on the function i
    redun.out <- getRedundancy(vars = function_names[i], species = species_names, data = data)
    sp.effect.aic[[i]] <- sapply(redun.out, function(x)(x))
    
  }
  
  # prepare the output
  sp.effect.aic <- as.data.frame(do.call(cbind, sp.effect.aic))
  names(sp.effect.aic) <- function_names
  sp.effect.aic <- cbind(species = rownames(sp.effect.aic), data.frame(sp.effect.aic, row.names = NULL))
  
  # return the data.frame with species effects
  return(sp.effect.aic)
  
}

# AIC_sp(data = mf.pa, function_names = f.names, species_names = spp.present)


# 2. calculate species effects using Gotelli et al. (2011)'s method

# data: data.frame with functions and species presence-absences as columns and samples as rows
# function_names: vector of function names
# species_names: vector of species names
# n_ran: number of randomisations

SES_score <- function(data, function_names, species_names, n_ran) {
  
  if(! "dplyr" %in% installed.packages()[,1]) stop(
    "this function requires the dplyr package to be installed"
  )
  if(! "tidyr" %in% installed.packages()[,1]) stop(
    "this function requires the tidyr package to be installed"
  )
  
  # load dplyr and tidyr packages
  library(dplyr)
  library(tidyr)
  
  # randomise the data and calculate d-values for each randomisation
  ran.d.vals <- vector("list", length = n_ran)
  for (i in 1:n_ran) {
    
    # randomise function values across plots
    data.random <- row.randomiser(func.names = function_names,
                                  species.names = species_names, 
                                  adf.data = data)
    
    # calculate difference values for each species for each function
    ran.d.vals[[i]] <-  
      data.random %>%
      pivot_longer(cols = all_of(species_names),
                   names_to = "species",
                   values_to = "pa") %>%
      group_by(species, pa) %>%
      summarise(across(.cols = all_of(function_names), mean ), .groups = "drop") %>%
      group_by(species) %>%
      summarise(across(.cols = all_of(function_names), diff ), .groups = "drop")
    
  }
  
  # bind the list into a data.frame  
  ran.d.vals <- bind_rows(ran.d.vals, .id = "run")
  
  # calculate mean and sd of d-values of randomised data
  ran.d.vals <- 
    ran.d.vals %>%
    pivot_longer(cols = all_of(function_names),
                 names_to = "function_name",
                 values_to = "d_value") %>%
    arrange(species, function_name) %>%
    group_by(species, function_name) %>%
    summarise(mean_d_value = mean(d_value, na.rm = TRUE),
              sd_d_value = sd(d_value, na.rm = TRUE), .groups = "drop")
  
  # calculate observed d values
  obs.d.vals <- 
    data %>%
    pivot_longer(cols = all_of(species_names),
                 names_to = "species",
                 values_to = "pa") %>%
    group_by(species, pa) %>%
    summarise(across(.cols = all_of(function_names), mean ), .groups = "drop") %>%
    group_by(species) %>%
    summarise(across(.cols = all_of(function_names), diff ), .groups = "drop") %>%
    pivot_longer(cols = starts_with("F_"),
                 names_to = "function_name",
                 values_to = "d_value_obs") %>%
    arrange(species, function_name)
  
  # join the observed d.vals to the randomised d.vals
  d.val.dat <- full_join(ran.d.vals, obs.d.vals, by = c("species", "function_name"))
  
  # calculate species importance scores
  sp.SES <- 
    d.val.dat %>%
    mutate(SES = (d_value_obs - mean_d_value)/(sd_d_value)) %>%
    mutate(SES_effect = if_else(SES > 2, 1, 
                                if_else(SES < -2, -1, 0))) %>%
    select(species, function_name, SES_effect) %>%
    pivot_wider(id_cols = "species",
                names_from = "function_name",
                values_from = "SES_effect")
  
  return(sp.SES)
  
}


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
        func.in <- function_matrix[function_matrix$species %in% sp.present, ]
        
        
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

### END
