
# Project: Review of multifunctionality in ecology, conservation and ecosystem service science

# Title: Does the number of functions matter?

# load relevant libraries
library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)
library(here)

# load functions from important scripts
source(here("Scripts/MF_functions_collated.R"))


# BEF slope and n-functions

# get the list of matrices of function combinations
function.combinations <- function(vector.func.names) {
  func.combs <- vector("list", length = (length(vector.func.names)-1) )
  for (i in 2:length(vector.func.names)){
    func.combs[[i-1]] <- combn(x = vector.func.names, m = i)
  }
  return(func.combs)
}

# write a function to flatten an individual part of a nested list
flatten.list.matrices <- function(nested.list){
  unnested.list <- split(nested.list, col(nested.list)) 
  names(unnested.list) <- NULL 
  return(unnested.list)
}

# combine function.combinations and flatten.list.matrices and loop over each n-func
get.function.combinations <- function(function.names){
  
  # get list of matrices with function combinations
  list.func.matrix <- function.combinations(vector.func.names = function.names)
  
  # flatten the first matrix in the list
  list.combination <- flatten.list.matrices(nested.list = list.func.matrix[[1]])
  
  # loop over this and bind into a list
  for (i in 2:length(list.func.matrix)){
    x <- flatten.list.matrices(nested.list = list.func.matrix[[i]])
    list.combination <- c(list.combination, x)
  }
  return(list.combination)
  
}


# function to standardise functions and translate them by minimum absolute value
standardise <- function(x) {
  mean.x <- mean(x)
  sd.x <- sd(x)
  x.standardised = ((x-mean.x)/sd.x)
  x.standardised.positive = x.standardised + abs(min(x.standardised))
  return(x.standardised.positive)
}


# using the defined functions:

# set up a function to calculate the realised diversity-multifunctionality slope
# for each combination of a set of ecosystem functions

# adf.data: data.frame where rows are plots and which contains ecosystem functions of interest
# mf.func.names: vector of function names to consider

get_BEF_mf_est <- function(adf.data, mf.func.names) {
  
  if (length(mf.func.names) < 2) {
    stop("must have more than two functions to adequately examine the effect of changing the number of functions")
  }
  
  # use the function to get a list of all function combinations
  list.func.names <- get.function.combinations(function.names = mf.func.names)
  
  # set an output list
  list.out <- vector("list", length(list.func.names))
  
  # loop over each set of function names
  for (i in 1:length(list.func.names)) {
    
    # get the vector of function names
    sample.func.names <- list.func.names[[i]]
    
    # subset out the functions and standardise them
    adf.func <- 
      adf.data %>%
      select(all_of(sample.func.names )) %>%
      mutate(across(.cols = all_of(sample.func.names) , ~standardise(.) ) )
    
    # calculate the multifunctionality metrics
    data.mf <- 
      adf.data %>%
      mutate(`sum MF` = MF_sum(adf = adf.func, vars = sample.func.names),
             `ave. MF` = MF_av(adf = adf.func, vars = sample.func.names),
             `SAM MF` = MF_dooley(adf = adf.func, vars = sample.func.names),
             `ENF MF` = as.numeric(hill_multifunc(adf = adf.func, vars = sample.func.names, scale = 1, HILL = TRUE)),
             `thresh.30 MF` = single_threshold_mf(adf = adf.func, vars = sample.func.names, thresh = 0.3),
             `thresh.70 MF` = single_threshold_mf(adf = adf.func, vars = sample.func.names, thresh = 0.7))
    
    # subset the multifunctionality metric names
    mf.names <- names(data.mf)[grepl(pattern = " MF", x = names(data.mf))]
    
    # for each multifunctionality metric, calculate the BEF-slope
    bef_mf_slope <- 
      sapply(mf.names, function(x){
        lm.x <- lm(reformulate("realised_diversity", as.name(x) ), data = data.mf)
        lm.x$coefficients[2]
      })
    names(bef_mf_slope) <- NULL # remove names from the vector
    
    # bind this into a data.frame
    df.out <- 
      data.frame(sim.id = i,
                 number_of_functions = length(sample.func.names),
                 multifunctionality_metric = mf.names,
                 realised_diversity_mf_est = bef_mf_slope)
    
    list.out[[i]] <- df.out
  }
  
  bind_rows(list.out)
  
}


### experiment on the Jena data

# load the cleaned Jena data
jena.dat <- read_csv(file = here("data/jena_data_cleaned.csv"))

# define variable groups
var.names <- names(jena.dat)

# (1) get species names
spp.p <- ( grepl("+\\.+", var.names) & nchar(var.names) == 7 )
spp.names <- var.names[spp.p]

# (2) get site identifiers
site.id <- c("year", "sowndiv", "plotcode", "realised_diversity")

# (3) get function names
func.names <- var.names[!(var.names %in% ( c(spp.names, site.id) ))  ]


# test this function
get_BEF_mf_est(adf.data = select(jena.dat, all_of(c(site.id, func.names))) , 
               mf.func.names = func.names[1:3])












