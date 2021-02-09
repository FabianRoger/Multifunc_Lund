
# Project: Review of multifunctionality in ecology, conservation and ecosystem service science

# Title: Implement the turnover approach using the multifunc package and an analogue to Meyer et al. (2017)'s null model

# install and load the multifunc package
# library(devtools)
# install_github("jebyrnes/multifunc")
library(multifunc)

# arguments for the functions:

# vars: vector of function names
# species: vector of species names
# data: data.frame with species abundances/presence-absences and function data
# output: either "prop_species" for proportion of species affecting a function
# or "overlap" which reports the mean pairwise overlap in species positive and negative effects
# between all functions

# n: number of null replicates

turnover_aic <- function(func.names, species.names, adf.data, output = "prop_species") {
  
  if(! "multifunc" %in% installed.packages()[,1]) stop(
    "this function requires the multifunc package to be installed"
  )
  
  # load the multifunc package
  library(multifunc)
  
  # use getRedundancy to get the effect of different species on each function
  redund.dat <- getRedundancy(vars = func.names, species = species.names, data = adf.data)
  
  # specify what to output
  if (output == "prop_species") {
    
    # calculate the number of species that positively affect different numbers of functions
    posCurve <- divNeeded(redund.dat, type = "positive")
    posCurve$div <- posCurve$div/ncol(redund.dat)
    row.names(posCurve) <- NULL
    posCurve$effect_direction <- "positive"
    
    # calculate the number of species that negative affect different numbers of functions
    negCurve<-divNeeded(redund.dat, type = "negative")
    negCurve$div<-negCurve$div/ncol(redund.dat)
    row.names(negCurve) <- NULL
    negCurve$effect_direction <- "negative"
    
    # bind these positive and negative effects together
    species_effects <- rbind(posCurve, negCurve)
    species_effects <- cbind(row.id = 1:nrow(species_effects), species_effects)
    
    return(species_effects)
  }
  
  else if(output == "overlap") {
    
    # calculate average sorensen overlap among all pairs of functions
    
    # positive function effects
    pos.overlap <- getOverlapSummary(redund.dat, 
                                     m=2, type = "positive", 
                                     index = "sorensen", denom = "set")[1]
    
    # negative function effects
    neg.overlap <- getOverlapSummary(redund.dat, m=2, 
                                     type = "negative", 
                                     index = "sorensen", denom = "set")[1]
    
    # bind this into a data.frame
    overlap_effects <- 
      data.frame(direction = c("positive", "negative"),
                 mean_sorensen_overlap = c(pos.overlap, neg.overlap))
    
    return(overlap_effects)
  } 
  
  else { print("error, specify an output") }
  
}

# function to mix-up relationship between function and species
row.randomiser <- function(func.names, species.names, adf.data) {
  
  # subset a matrix of functions
  func.mat <- adf.data[, func.names]
  func.mat.nrow <- nrow(func.mat) # calculat the number of rows in that matrix
  
  # use sample to get random row ids from the function matrix
  random.row.ids <- sample(x = 1:func.mat.nrow , size = func.mat.nrow , replace = FALSE)
  
  # randomise the function matrix
  func.mat.random <- func.mat[random.row.ids, ]
  
  # get a matrix of species data
  spec.mat <- adf.data[, species.names]
  
  # bind this randomised data together
  adf.data.random <- cbind(spec.mat, func.mat.random)
  
  return(adf.data.random)
  
}


# use the defined functions to calculate turnover for randomised data and observed data
turnover_aic_null <- function(func.names, 
                              species.names, 
                              adf.data, 
                              output = "prop_species",
                              n = 100) {
  
  if(! "dplyr" %in% installed.packages()[,1]) stop(
    "this function requires the dplyr package to be installed"
  )
  
  # randomise the data and calculate turnover indices n times
  null.out <- vector("list", length = n)
  for (i in 1:n){
    
    # randomise the function positions relative to the species using row.randomiser
    data.random <- row.randomiser(func.names = func.names, 
                                  species.names = species.names, 
                                  adf.data = adf.data)
    
    # apply the approach to the randomised data
    turnover.out <- turnover_aic(func.names = func.names, 
                                 species.names = species.names, 
                                 adf.data = data.random,
                                 output = output)
    
    # write the n_func proportion of the species pool to a list
    null.out[[i]] <- turnover.out
  }
  
  # bind the loop output into a data.frame
  null.out <- dplyr::bind_rows(null.out, .id = "null.rep")
  
  # get observed value
  obs.out <- turnover_aic(func.names = func.names, 
                          species.names = species.names, 
                          adf.data = data.random,
                          output = output)
  
  # bind this into a list
  return(list(null.out, obs.out))
  
}


# test the turnover_aic_null function and plot the results
x <- turnover_aic_null(func.names = vars, 
                              species.names = species, 
                              adf.data = germany,
                              output = "prop_species",
                              n = 1000)


# plot the null expectations versus the observed data
library(ggplot2)
ggplot() +
  stat_smooth(data = x[[1]],
              mapping = aes(x = nfunc, y = div, group = null.rep),
              geom='line', alpha=0.01, size = 1, se=FALSE, method = "lm") +
  geom_smooth(data = x[[2]],
              mapping = aes(x = nfunc, y = div), colour = "red",
              method = "lm", se = FALSE) +
  facet_wrap(~effect_direction, scales = "free") +
  theme_classic()



























