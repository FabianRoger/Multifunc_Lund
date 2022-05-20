
#' @title theme_meta
#' @description Function to efficiently output the slope of a linear model
#' 
#' @author James G. Hagan (james_hagan(at)outlook.com)
#' @param data data.frame with data to be used in the linear model
#' @param response string containing the name of the response variable (y-variable)
#' @param explanatory string containing the name of the explanatory variable (x-variable)
#' 
#' @return returns a data.frame with a single named column ("explanatory_response") and the slope estimate
#' 


lm.cleaner <- function(data, response, explanatory) {
  
  # check that the correct packages are installed
  if(! "dplyr" %in% installed.packages()[,1]) stop(
    "this function requires dplyr to be installed"
  )
  
  # check that the correct packages are installed
  if(! "broom" %in% installed.packages()[,1]) stop(
    "this function requires broom to be installed"
  )
  
  # load the dplyr and broom libraries
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
