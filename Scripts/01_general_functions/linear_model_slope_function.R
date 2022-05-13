
# Project: Review of multifunctionality in ecology, conservation and ecosystem service science

# Title: Function to efficiently output the slope of a linear model

# arguments

# data - data.frame with data to be used in the linear model
# response - name of the response variable
# explanatory - name of the explanatory variable

# define function to efficiently output the slope
lm.cleaner <- function(data, response, explanatory) {
  
  # check that the correct packages are installed
  if(! "dplyr" %in% installed.packages()[,1]) stop(
    "this function requires dplyr to be installed"
  )
  
  # check that the correct packages are installed
  if(! "broom" %in% installed.packages()[,1]) stop(
    "this function requires broom to be installed"
  )
  
  library(dplyr)
  library(broom)
  
  lm.x <- 
    lm(reformulate(explanatory, response), data = data) %>% 
    tidy %>% 
    filter(term == explanatory) %>% 
    dplyr::select(estimate)
  
  names(lm.x) <- paste(explanatory, response, sep = "_")
  
  return(lm.x)
  
}

### END
