
# Project: Multifunctionality workshop

# Title: PCA approach analysis

library(ggplot2)
library(dplyr)
library(tidyr)
library(here)

source(here("Scripts/Multifunctionality-Simulations/Multifunc_simulations_functions.R"))
source(here("Scripts/MF_functions_collated.R"))
source(here("Scripts/function_plotting_theme.R"))


# set seed to replicate analysis: set.seed(1600436230)


# set number of runs to do
n <- 10

# set up a blank list to fill
pca_func <- vector("list", length = n)

# run the loop to generate three datasets with different correlation structures
for (i in (1:n) ) {
  
  set.seed(1555 + i)
  
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
  Sigma[1:2,1:2] <- 0.7
  Sigma[4:6,4:6] <- 0.5
  Sigma[7:9,7:9] <- 0.7
  
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
  pca_func[[i]] <- AvFunc
  
}


# calculate the MF metrics on these datasets with different levels of correlation

pca_mf <- 
  
  lapply(pca_func, function(x) {
    
    mutate(x, 
           Meyer_mf = pca_multifunc(adf = x, vars = func.names, standardise = FALSE),
           Av_mf = MF_av(adf = x, vars = func.names) )
    
  }
  
  )

pca_mf <- bind_rows(pca_mf, .id = "run")

ggplot(data = pca_mf,
       mapping = aes(x = Av_mf, y = Meyer_mf, colour = run)) +
  geom_point(alpha = 0.3, shape = 16) +
  geom_smooth(method = "lm", se = FALSE, size = 0.75) +
  scale_colour_viridis_d(option = "C", end = 0.9) +
  theme_meta() +
  theme(legend.position = "none")


