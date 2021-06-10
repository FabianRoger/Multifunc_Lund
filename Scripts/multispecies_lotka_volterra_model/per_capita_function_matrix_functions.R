
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
  func.mat <- cbind(species_list, func.mat)
  names(func.mat) <- c("species", func.names)
  
  return(func.mat)
  
}


# example:

# func_matrix_generator(species_list = c("sp_1", "sp_2"), 
                      # func.n = 5, func.spec = "specialist", 
                      # prob.neg = 0)

### END
