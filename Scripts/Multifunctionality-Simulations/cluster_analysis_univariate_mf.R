
# Project: Multifunctionality workshop

# Title: Cluster analysis of different univariate multifunctionality metrics

library(ggplot2)
library(dplyr)
library(tidyr)
library(corrplot)
library(here)

source(here("Scripts/Multifunctionality-Simulations/Multifunc_simulations_functions.R"))
source(here("Scripts/MF_functions_collated.R"))


# set seed to replicate analysis: set.seed(1600436230)

# set up correlations among the three function clusters
cor_struc <- 
  list(c1 = c(0, 0, 0),
       c2 = c(0.7, 0.5, 0),
       c3 = c(0.7, 0.5, 0.7))

# set up a blank list to fill
clu_func <- vector("list", length = length(cor_struc))

# run the loop to generate three datasets with different correlation structures
for (i in 1:length(clu_func) ) {
  
  set.seed(1600436230 + i)
  
  # number of species
  specnum <- 10
  
  # number of functions
  funcnum <- 9
  
  # distribution from which to draw function values
  # argument is irrelevant here but needs to be set
  distribution = "runif"
  
  # choose pairwise correlation strength
  COR <- 0
  
  # make correlation matrix (strictly speaking a covariance matrix but for these simulations it does not matter)
  Sigma <- matrix(COR, ncol = funcnum, nrow = funcnum)
  
  # make three 'cluster' of correlated functions
  Sigma[1:2,1:2] <- cor_struc[[i]][1]
  Sigma[4:6,4:6] <- cor_struc[[i]][2]
  Sigma[7:9,7:9] <- cor_struc[[i]][3]
  
  diag(Sigma) <- 1
  
  Sigma
  
  # create function matrix
  FuncMat <- FunctionValue(specnum,funcnum, distribution, min = 0, max = 1)
  FuncMat
  
  # replace function values with correlated function values
  FuncMat_long <- 
    FuncMat %>% 
    pivot_wider(names_from = Functions, 
                values_from = Funcval)
  
  # draw correlated functions (with mean 0)
  corF <- mvrnorm(n = specnum, mu = rep(0, funcnum), Sigma = Sigma)
  
  # shift to positive
  corF <- apply(corF, 2, function(x){ x + abs(min(x)) })
  
  # backtransform to long dataset
  FuncMat_long[, 2:(funcnum+1) ] <- corF
  
  FuncMat <-
    FuncMat_long %>% 
    pivot_longer(cols = starts_with("F"),
                 names_to = "Functions",
                 values_to = "Funcval") 
  
  # extract function names
  func.names <- as.character( unique( FuncMat$Functions))
  
  
  # use these function values to generate data from a biodiversity experiment
  maxrep <- choose(specnum, floor(specnum/2))
  
  # simulate plot x species matrix
  SpecMat <- SpeciesMatrix(specnum = specnum, maxrep = maxrep)
  
  # select method and additional parameters if appropriate by setting the `if` statement to `TRUE`
  if (TRUE) {method = "av"}
  
  if (FALSE) {method = "comp"
  CF = 3
  compfunc = c("F 1", "F 6")
  r = 0.25}
  
  # Average function
  AvFunc <- AverageFunction(SpecMat, FuncMat,
                            method = method, 
                            CF = CF, 
                            compfunc = compfunc,
                            r = r)
  
  # add small normal error
  AvFunc <- 
    AvFunc %>% 
    mutate_at(func.names, function(x) {x + runif(n(), 0, 0.05)})
  
  # standardize by z-score
  AvFunc[,func.names] <- 
    
    apply(AvFunc[,func.names], 2, function(x) {
    
      y <- (x - mean(x))/sd(x)
      
      return(y+abs(min(y)))
  }
  
  )
  
  # add this output into an object
  clu_func[[i]] <- AvFunc
  
}


# calculate the MF metrics on these datasets with different levels of correlation

AvFunc_MF <- 
  
  lapply(clu_func, function(x) {
  
    mutate(x, 
           Hill_mf = hill_multifunc(adf = x, vars = func.names, scale = 1, HILL = TRUE),
           Manning30_mf = manning_multifunc(adf = x, vars = func.names, thresh = 0.3),
           Manning50_mf = manning_multifunc(adf = x, vars = func.names, thresh = 0.5),
           Manning70_mf = manning_multifunc(adf = x, vars = func.names, thresh = 0.7),
           Meyer_mf = pca_multifunc(adf = x, vars = func.names, standardise = FALSE),
           Pasari_mf = MF_pasari(adf = x, vars = func.names),
           Dooley_mf = MF_dooley(adf = x, vars = func.names),
           Jing_mf = MF_jing(adf = x, vars = func.names),
           Sum_mf = MF_sum(adf = x, vars = func.names),
           Av_mf = MF_av(adf = x, vars = func.names),
           Thresh30_mf = single_threshold_mf(adf = x, vars = func.names, thresh = 0.3),
           Thresh50_mf = single_threshold_mf(adf = x, vars = func.names, thresh = 0.5),
           Thresh70_mf = single_threshold_mf(adf = x, vars = func.names, thresh = 0.7),
           Mesli_mf = MF_mesli(adf = x, vars = func.names),
           Simpson_mf = MF_simpsons_div(adf = x, vars = func.names) )
  
}

)

AvFunc_MF








