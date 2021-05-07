
# Project: Review of multifunctionality in ecology, conservation and ecosystem service science

# Title: Code to generate per capita function coefficients for a list of species

# Relevant function is: func_matrix_generator

# parameter definitions

# species_list: vector of species names
# func.n: number of functions to simulate
# func.spec: simulates functional specialisation or generalisation
# - "specialist": gamma distribution (shape = 1, rate = 2)
# - "generalist": normal distribution (mean = 0.5, sd = 0.1)
# prob.neg: proportion of negative coefficients


# function to get functional coefficients
per_capita_functions <- function(func.n = 5, func.spec = "specialist", prob.neg = 0.1) {
  
  if(func.spec == "specialist") {
    
    f.coef <- rgamma(n = func.n, shape = 1, rate = 2)
    
  } else if (func.spec == "generalist") {
    
    f.coef <- rnorm(n = func.n, mean = 0.5, sd = 0.1)
    
  } else { stop("error! specify relevant specialist or generalist") }
  
  x <- f.coef - quantile(f.coef, prob.neg)
  y <- round(x, digits = 4)
  return(y)
  
}

# function to get function matrix for each species
func_matrix_generator <- 
  function(species_list, func.n, func.spec, prob.neg) {
  
  func.mat <- 
    lapply(species_list, function(x) {
      
      per_capita_functions(func.n = func.n, 
                           func.spec = func.spec, 
                           prob.neg = prob.neg)
      
    })
  
  # make a function name list
  func.names <- paste("F_", 1:func.n, sep = "")
  
  # bind this into a data.frame
  func.mat <- data.frame(do.call(rbind, func.mat))
  func.mat <- cbind(spp.list, func.mat)
  names(func.mat) <- c("species", func.names)
  
  return(func.mat)
  
}

# example:
func_matrix_generator(species_list = species, 
                      func.n = 5, func.spec = "specialist", prob.neg = 0)




### potentially tailor these distributions to empirical data...

### fitting distributions for the per capita functions

# load multifunc package
library(multifunc)

data(all_biodepth)
all_biodepth %>% View()
allVars<-qw(biomassY3, root3, N.g.m2,  light3, N.Soil, wood3, cotton3)
varIdx<-which(names(all_biodepth) %in% allVars)

biodepth_list <- split(all_biodepth, all_biodepth$location)

germany <- biodepth_list[[1]]

germany<-subset(all_biodepth, all_biodepth$location=="Germany")
vars<-whichVars(germany, allVars)
species<-relevantSp(germany,26:ncol(germany))
spIDX <- which(names(germany) %in% species) 

germany

mono <- subset(germany, rowSums(germany[,spIDX])==1)
names(mono)

mono <- 
  mono %>%
  dplyr::select(plot, all_of(vars))

library(tidyr)
mono_pc <- 
  mono %>%
  pivot_longer(cols = c("root3", "N.g.m2", "N.Soil", "cotton3"),
               names_to = "Func",
               values_to = "value") %>%
  mutate(per_capita_function = value/biomassY3) %>%
  dplyr::select(-value) %>%
  pivot_wider(names_from = "Func",
              values_from = "per_capita_function")

install.packages("fitdistrplus")
library(fitdistrplus)

# need to choose better distributions
dist_choose <- function(data, dist_list = c("weibull", "norm", "gamma", "cauchy", "lnorm")) {
  
  # load the fitdistrplus package
  library(fitdistrplus)
  
  aic_out <- vector(length = length(dist_list))
  for (i in 1:length(dist_list)) {
    
    x <- fitdist(data, distr = dist_list[i])
    aic_out[i] <- x$aic
    
  }
  
  # add names to the distribution
  names(aic_out) <- dist_list
  
  # select the lowest AIC distribution
  low_aic_dist <- names(sort(aic_out, decreasing = FALSE)[1])
  
  y <- fitdist(data, distr = low_aic_dist)
  df <- data.frame(y$estimate)
  row.names(df)
  
  df$parameter <- row.names(df)
  
  row.names(df) <- NULL
  
  df$distribution <- low_aic_dist
  
  df <- 
    df %>%
    dplyr::select(distribution, parameter, y.estimate)
  
  return(df)
  
}

x <- dist_choose(data = mono_pc$cotton3, 
            dist_list = c("gamma"))
x

a <- rgamma(n = 1000, shape = 1, rate = 2)
hist(a)
mean(a)

b <- rnorm(n = 1000, mean = 0.5, sd = 0.1)

hist(b)
mean(b)

dist_choose(data = mono_pc$N.g.m2, 
            dist_list = c("weibull"))

dist_choose(data = mono_pc$N.g.m2, 
            dist_list = c("cauchy"))
