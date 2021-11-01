
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
  return(as_tibble(sp.effect.aic))
  
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
    pivot_longer(cols = all_of(function_names),
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

# function to calculate proportion of the species pool required to drive function

# data: data.frame with functions and species presence-absences as columns and samples as rows
# function_names: vector of function names
# species_names: vector of species names
# method: method to calculate species effects ("AIC" or "SES")
# n_ran: number of randomisations

prop_species_pool <- function(data, function_names, species_names, method = "AIC", n_ran = 100) {
  
  # get the function combinations
  f.combs <- c(function_names, get.function.combinations(function.names = function_names))
  
  # set up the list
  prop_sp_out <- vector("list", length = length(f.combs))
  
  # for each function combination, calculate the proportion of the species pool 
  # that has at least one positive effect using either:
  # aic-based turnover approach
  # SES-based turnover approach
  
  for(i in 1:length(f.combs)) {
    
    if (method == "AIC") {
      
      df.in <- AIC_sp(data = data, function_names = f.combs[[i]], species_names = species_names)
      
    } else if (method == "SES") {
      
      df.in <- SES_score(data = data, function_names = f.combs[[i]], species_names = species_names, n_ran = n_ran)
      
    } else { stop("choose appropriate method for calculating the species pool") }
    
    # test if each species contributes positively to at least one function
    x <- apply(df.in[ , -1], 1, function(z) { ifelse(any(z > 0), TRUE, FALSE) })
    
    # count the number of species positively contributing to each function and divide by total number of species
    prop_pos <- length(df.in[x, ]$species)/length(species_names)
    
    prop_sp_out[[i]] <- tibble(number_functions = length(f.combs[[i]]),
                               function_id = paste(f.combs[[i]], collapse = "."),
                               prop_species_pool = prop_pos)
    
  }
  
  dplyr::bind_rows(prop_sp_out)
  
}

# create a function to randomise the data and then calculate the proportion
# of species pool

# data: data.frame with functions and species presence-absences as columns and samples as rows
# function_names: vector of function names
# species_names: vector of species names
# method: method to calculate species effects ("AIC" or "SES")
# n_ran: number of randomisations

prop_species_pool_random <- function(data, function_names, species_names, method = "AIC", n_ran = 100, n = 10) {

  # create n random datasets in a list
  random_rows <- vector("list", length = n)
  
  for (i in 1:n) {
    
    random_rows[[i]] <- row.randomiser(func.names = function_names,
                                       species.names = species_names,
                                       adf.data = data)
    
  }
  
  # apply over this list
  prop_list<- 
    lapply(random_rows, function(x) {
      
      prop_species_pool(data = x, function_names = function_names,
                        species_names = species_names, method = method, n_ran = n_ran)
      
    })
  
  # bind this into a data.frame
  bind_rows(prop_list, .id = "run")
  
  
}

### END
