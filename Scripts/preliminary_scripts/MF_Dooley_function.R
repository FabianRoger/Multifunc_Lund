
# function to calculate three multifunctionality metric proposed in the literature:
# Dooley (2018) Scaled Average Multifunctionality (SAM)

# adf, is dataframe with plots in rows, and functions in columns
# vars has to bee a named vector of functions to include which has to correspond to column names
# stand is the method of standardisation. three methods are supported: (1) no standardisation ("none") (2) z-score standardisation ("z_score") and (3) by the maximum ("max")

MF_dooley <- function(adf, vars, stand = "none") {
  
  # extract functions from input matrix: adf
  adf_mat <- adf[, vars] 
  
  # standardise the functions using: (1) the z-score ("z_score") or (2) by the maximum ("max")
  
  if (stand == "none") {
    
  } else if (stand == "z_score") {
    
    adf_mat <- as.matrix(scale(x = adf_mat, center = TRUE, scale = TRUE))
    
  } else if (stand == "max") {
    
    adf_mat <- 
      apply(X = adf_mat, MARGIN = 2, FUN = function(z) { y <- z/max(z)} )
    
  } else { 
    
    stop("this function requires a standardisation method") 
    
  }
  
  adf$mf_dooley <- apply(adf_mat, MARGIN = 1, function(x) mean(x) )/apply(adf_mat, MARGIN = 1, function(x) sd(x) )
  
  return(adf)
  
  }
  
  
  
  
  