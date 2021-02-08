
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


# BEF slope and n-functions

# get the list of matrices of function combinations
function.combinations <- function(vector.func.names) {
  func.combs <- vector("list", length = length(func.names))
  for (i in n.funcs){
    func.combs[[i]] <- combn(x = vector.func.names, m = i)
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

# use the function to get a list of all function combinations
get.function.combinations(function.names = func.names)














