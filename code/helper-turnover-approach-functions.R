#'
#' @title Functions to implement turnover approaches
#' 
#' @description This script contains a series of functions that are used
#' to implement the turnover (aka overlap) approach to multifunctionality using
#' two separate methods: Hector and Bagchi's (2007) information theoretic
#' approach (implemented as per Byrnes et al. 2014) and a novel implementation
#' that uses Gotelli et al.'s (2011) species importance scores instead of the
#' information theoretic approach.
#'

#' @title row_randomiser
#' 
#' @param func_names - vector of function names
#' @param sp_names - vector of species names
#' @param adf_data - dataset containing functions and species columns
#' 

# function to mix-up relationship between function and species
row_randomiser <- function(func_names, sp_names, adf_data) {
  
  assertthat::assert_that(
    is.vector(func_names) && is.vector(sp_names),
    msg = paste("either func_names or sp_names are not vectors")
  )
  
  assertthat::assert_that(
    is.data.frame(adf_data) | dplyr::is.tbl(adf_data),
    msg = paste("adf_data is not a data.frame or tibble")
  )
  
  assertthat::assert_that(
    all(func_names %in% names(adf_data)) && all(sp_names %in% names(adf_data)),
    msg = "not all function names nor all species names are present in the input data"
  )
  
  assertthat::assert_that(
    length(func_names) > 1,
    msg = paste("func_names must have more than one function")
  )
  
  assertthat::assert_that(
    length(sp_names) > 1,
    msg = paste("sp_names must have more than one species")
  )
  
  # subset a matrix of functions
  func_mat <- adf_data[, func_names]
  
  # calculate the number of rows in that matrix
  func_mat_nrow <- nrow(func_mat)
  
  # use sample to get random row ids from the function matrix
  random_row_ids <- sample(x = 1:func_mat_nrow, size = func_mat_nrow , replace = FALSE)
  
  # randomise the function matrix
  func_mat_random <- func_mat[random_row_ids, ]
  
  # get a matrix of species data
  spec_mat <- adf_data[, sp_names]
  
  # bind this randomised data together
  adf_data_random <- cbind(spec_mat, func_mat_random)
  
  return( as.data.frame(adf_data_random) )
  
}

#' @title function_combinations
#' 
#' @param func_names - vector of function names
#'

# get the list of matrices of function combinations
function_combinations <- function(func_names) {
  
  assertthat::assert_that(
    is.vector(func_names),
    msg = paste("func_names is not a vector")
  )
  
  assertthat::assert_that(
    length(func_names) > 1,
    msg = paste("func_names must have more than one function")
  )
  
  nested_list_matrices <- vector("list", length = (length(func_names)-1) )
  for (i in 2:length(func_names)){
    nested_list_matrices[[i-1]] <- combn(x = func_names, m = i)
  }
  return(nested_list_matrices)
  
}

#' @title flatten_list_matrices
#' 
#' @param nested_list_matrices - output from the function_combinations() function
#'
#' @description the function_combinations() function outputs a nested list of
#' function combinations i.e. combinations of two functions is one nested layer
#' and combinations of three functions is another nested layer. This function
#' converts those nested layers into a single list where each list element
#' is a different function combination
#'

flatten_list_matrices <- function(nested_list_matrices){
  
  unnested_list <- split(nested_list_matrices, col(nested_list_matrices)) 
  names(unnested_list) <- NULL 
  return(unnested_list)
  
}

#' @title get_function_combinations
#' 
#' @param func_names - vector of function names
#'
#' @description function combines several functions to output a list where
#' each element is a combination of functions.
#'

get_function_combinations <- function(func_names){
  
  # get raw list of matrices
  list_func_matrix <- function_combinations(func_names = func_names)
  
  if (length(func_names) < 3) {
    
    # flatten the first matrix in the list
    list_combination <- flatten_list_matrices(nested_list_matrices = list_func_matrix[[1]])
    
  } else {
    
    # flatten the first matrix in the list
    list_combination <- flatten_list_matrices(nested_list_matrices = list_func_matrix[[1]])
    
    # loop over this and bind into a list
    for (i in 2:length(list_func_matrix)){
      x <- flatten_list_matrices(nested_list_matrices = list_func_matrix[[i]])
      list_combination <- c(list_combination, x)
    }
    
  }
  
  return(list_combination)
  
}

#' @title AIC_sp
#' 
#' @param data - data.frame containing functions and species abundances or PAs
#' @param func_names - vector of function names
#' @param sp_names - vector of species names
#'
#' @description use the AIC-approach as implemented by Byrnes et al. (2014)
#' in the multifunc package to estimate the number of species required to 
#' support a set of functions
#'

AIC_sp <- function(data, func_names, sp_names) {
  
  assertthat::assert_that(
    is.vector(func_names) && is.vector(sp_names),
    msg = paste("either func_names or sp_names are not vectors")
  )
  
  assertthat::assert_that(
    is.data.frame(data) | dplyr::is.tbl(data),
    msg = paste("data is not a data.frame or tibble")
  )
  
  assertthat::assert_that(
    all(func_names %in% names(data)) && all(sp_names %in% names(data)),
    msg = "not all function names nor all species names are present in the input data"
  )
  
  # make sure the multifunc package is installed
  if(! "multifunc" %in% installed.packages()[,1]) stop(
    "this function requires the multifunc package to be installed"
  )
  
  # make sure the dplyr package is installed
  if(! "dplyr" %in% installed.packages()[,1]) stop(
    "this function requires the dplyr package to be installed"
  )
  
  sp_effect_aic <- vector("list", length = length(func_names))
  for (i in 1:length(func_names)){
    
    # get species effect on each function using abundances
    # this outputs a vector of species effects (-1, 0 or 1) on the function i
    redun_out <- multifunc::getRedundancy(vars = func_names[i], species = sp_names, data = data)
    sp_effect_aic[[i]] <- sapply(redun_out, function(x)(x))
    
  }
  
  # prepare the output
  sp_effect_aic <- as.data.frame(do.call(cbind, sp_effect_aic))
  names(sp_effect_aic) <- func_names
  sp_effect_aic <- cbind(species = rownames(sp_effect_aic), data.frame(sp_effect_aic, row.names = NULL))
  
  # return the data.frame with species effects
  return( dplyr::as_tibble(sp_effect_aic) )
  
}

#' @title SES_sp
#' 
#' @param data - data.frame containing functions and species abundances or PAs
#' @param func_names - vector of function names
#' @param sp_names - vector of species names
#' @param n_ran - number of randomisations
#'
#' @description use the SES-approach as implemented by Gotelli et al. (2011)
#' to estimate the number of species required to support a set of functions
#'

SES_sp <- function(data, func_names, sp_names, n_ran = 100) {
  
  assertthat::assert_that(
    is.vector(func_names) && is.vector(sp_names),
    msg = paste("either func_names or sp_names are not vectors")
  )
  
  assertthat::assert_that(
    is.data.frame(data) | dplyr::is.tbl(data),
    msg = paste("data is not a data.frame or tibble")
  )
  
  assertthat::assert_that(
    all(func_names %in% names(data)) && all(sp_names %in% names(data)),
    msg = "not all function names nor all species names are present in the input data"
  )
  
  assertthat::is.count(n_ran)
  
  sp_dat <- data[, sp_names]
  func_dat <- data[, func_names]
  
  func_list <- vector("list", length = length(func_names))
  for (i in 1:length(func_names)) {
    
    # get the focal function
    x <- func_dat[[func_names[i]]]
    
    sp_out <- 
      lapply(sp_dat, function(y) {
        
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
    
    func_list[[i]] <- unlist(sp_out)
    
  }
  
  # sort out the output
  sp_effect_ses <- do.call("cbind", func_list)
  sp_effect_ses <- as.data.frame(sp_effect_ses)
  names(sp_effect_ses) <- func_names
  sp_effect_ses <- dplyr::as_tibble(sp_effect_ses, rownames = "species")
  
  return(sp_effect_ses)
  
}

#' @title prop_species_pool
#' 
#' @param data - data.frame containing functions and species abundances or PAs
#' @param func_names - vector of function names
#' @param sp_names - vector of species names
#' @param method - AIC or SES to implement the different types of turnover approaches
#' @param n_ran - number of randomisations if SES is chosen
#'
#' @description calculates the proportion of
#'

prop_species_pool <- function(data, func_names, sp_names, method = "AIC", n_ran = 100) {
  
  # start by getting the effect of each species on the different functions using either:
  # 1. AIC-approach based on the multifunc package
  # 2. SES-approach of Gotelli et al. (2011)
  
  if (method == "AIC") {
    
    df.in <- AIC_sp(data = data, function_names = function_names, sp_names = sp_names)
    
  } else if (method == "SES") {
    
    df.in <- SES_score(data = data, function_names = function_names, sp_names = sp_names, n_ran = n_ran)
    
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
      
      prop_pos <- length(df.in[x, ]$species)/length(sp_names)
      
    }
    
    # tally up negative effects
    
    # test if each contributes negatively to at least one function
    x2 <- apply(y, 1, function(z) { ifelse(any(z < 0), TRUE, FALSE) })
    
    # if there are no species contributing to function then set prop_pos FALSE
    if (all(x2 == FALSE)) {
      
      prop_neg <- 0
      
      # count the number of species positively contributing to each function and divide by total number of species    
    } else {
      
      prop_neg <- length(df.in[x2, ]$species)/length(sp_names)
      
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
# sp_names: vector of species names
# method: method to calculate species effects ("AIC" or "SES")
# n_ran: number of randomisations when using the SES method
# n: number of random datasets

prop_species_pool_random <- function(data, function_names, sp_names, method = "AIC", n_ran = 100, n = 10) {
  
  # create n random datasets in a list
  random_rows <- vector("list", length = n)
  
  for (i in 1:n) {
    
    random_rows[[i]] <- row.randomiser(func.names = function_names,
                                       species.names = sp_names,
                                       adf.data = data)
    
  }
  
  # apply over this list
  prop_list<- 
    lapply(random_rows, function(x) {
      
      prop_species_pool(data = x, function_names = function_names,
                        sp_names = sp_names, method = method, n_ran = n_ran)
      
    })
  
  # bind this into a data.frame
  bind_rows(prop_list, .id = "run")
  
}

### END
