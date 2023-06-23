
# Project: Review of multifunctionality in ecology, conservation and ecosystem service science

# Title: Implement the turnover approach using the multifunc package and an analogue to Meyer et al. (2017)'s null model

# set-up general functions

# function to mix-up relationship between function and species
row.randomiser <- function(func.names, species.names, adf.data) {
  
  # subset a matrix of functions
  func.mat <- adf.data[, func.names]
  func.mat.nrow <- nrow(func.mat) # calculate the number of rows in that matrix
  
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

SES_score <- function(data, function_names, species_names, n_ran = 100) {
  
  sp.dat <- data[, species_names]
  func.dat <- data[, function_names]
  
  func.list <- vector("list", length = length(function_names))
  for (i in 1:length(function_names)) {
    
    x <- func.dat[[function_names[i]]]
    
    sp.out <- 
      lapply(sp.dat, function(y) {
        
        obs_D <- mean(x[y == 1]) - mean(x[y == 0])
        
        ran_D <- vector("double", length = n_ran)
        for (i in 1:n_ran) {
          
          z <- sample(x, size = length(x), replace = TRUE)
          ran_D[i] <- mean(z[y == 1]) - mean(z[y == 0])
          
        }
        
        SES_i <- (obs_D - mean(ran_D) )/sd(ran_D)
        
        # assign positive or negative effects depending on the effect size (>2 or <-2)
        if (SES_i > 2) {
          
          SES_i <- 1
          
        } else if (SES_i < -2) {
          
          SES_i <- -1 
          
        } 
        
        else { SES_i <- 0 }
        
      }
      
      )
    
    func.list[[i]] <- unlist(sp.out)
    
  }
  
  # sort out the output
  df <- do.call("cbind", func.list)
  df <- as.data.frame(df)
  names(df) <- function_names
  df <- as_tibble(df, rownames = "species")
  
  return(df)
  
}

# function to calculate proportion of the species pool required to drive function

# data: data.frame with functions and species presence-absences as columns and samples as rows
# function_names: vector of function names
# species_names: vector of species names
# method: method to calculate species effects ("AIC" or "SES")
# n_ran: number of randomisations when using the SES score

prop_species_pool <- function(data, function_names, species_names, method = "AIC", n_ran = 100) {
  
  # start by getting the effect of each species on the different functions using either:
  # 1. AIC-approach based on the multifunc package
  # 2. SES-approach of Gotelli et al. (2011)
  
  if (method == "AIC") {
    
    df.in <- AIC_sp(data = data, function_names = function_names, species_names = species_names)
    
  } else if (method == "SES") {
    
    df.in <- SES_score(data = data, function_names = function_names, species_names = species_names, n_ran = n_ran)
    
  } else { stop("choose appropriate method for calculating the species pool") }
  
  # get the function combinations
  f.combs <- c(function_names, get.function.combinations(function.names = function_names))
  
  # set up the list
  prop_sp_out <- vector("list", length = length(f.combs))
  
  # for each function combination, calculate the proportion of the species pool 
  # that has at least one positive effect
  
  for(i in 1:length(f.combs)) {
    
    # subset the functions of interest
    y <- df.in[, f.combs[[i]] ]
    
    # tally up positive effects
    
    # test if each species contributes positively to at least one function
    x <- apply(y, 1, function(z) { ifelse(any(z > 0), TRUE, FALSE) })
    
    # if there are no species contributing to function then set prop_pos FALSE
    if (all(x == FALSE)) {
      
      prop_pos <- 0
    
    # count the number of species positively contributing to each function and divide by total number of species    
    } else {
      
      prop_pos <- length(df.in[x, ]$species)/length(species_names)
      
    }
    
    # tally up negative effects
    
    # test if each contributes negatively to at least one function
    x2 <- apply(y, 1, function(z) { ifelse(any(z < 0), TRUE, FALSE) })
    
    # if there are no species contributing to function then set prop_pos FALSE
    if (all(x2 == FALSE)) {
      
      prop_neg <- 0
      
      # count the number of species positively contributing to each function and divide by total number of species    
    } else {
      
      prop_neg <- length(df.in[x2, ]$species)/length(species_names)
      
    }
    
    
    prop_sp_out[[i]] <- tibble(number_functions = length(f.combs[[i]]),
                               function_id = paste(f.combs[[i]], collapse = "."),
                               positive_effect = prop_pos,
                               negative_effect = prop_neg)
    
  }
  
  df.out <- dplyr::bind_rows(prop_sp_out)
  
  df.out <- tidyr::pivot_longer(data = df.out,
                                cols = c("positive_effect", "negative_effect"),
                                names_to = "effect_direction",
                                values_to = "proportion_species_pool")
  
  df.out <- dplyr::arrange(df.out, 
                           effect_direction, number_functions, function_id)
  
}

# create a function to randomise the data and then calculate the proportion
# of species pool

# data: data.frame with functions and species presence-absences as columns and samples as rows
# function_names: vector of function names
# species_names: vector of species names
# method: method to calculate species effects ("AIC" or "SES")
# n_ran: number of randomisations when using the SES method
# n: number of random datasets

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
