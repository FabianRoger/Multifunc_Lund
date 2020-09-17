
# function to calculate three multifunctionality metric proposed in the literature:
# (1) Pasari et al. (2013) multifunctionality metric
# (2) Dooley (2018) Scaled Average Multifunctionality (SAM)
# (3) Jing et al. (2020) Scaling Mulifunctionality Metric

# adf, is dataframe with plots in rows, and functions in columns
# vars has to bee a named vector of functions to include which has to correspond to column names
# stand is the method of standardisation. Two methods are supported: (1) z-score standardisation ("z_score") and (2) by the maximum ("max")
# metric is the metric to be calculated: there are (3) options to choose from:
# (1) Jing et al. (2020) : "scaling_mf"
# (2) Pasari et al. (2013): "unique_mf"
# (3) Dooley (2018): "sam_mf"

MF_scale_funcs <- function(adf, vars, 
                           stand = "z_score", 
                           metric = "scaling_mf") {
  
  # extract functions from input matrix: adf
  adf_mat <- adf[, vars] 
  
  # standardise the functions using: (1) the z-score ("z_score") or (2) by the maximum ("max")
  
  if (stand == "z_score") {
    
    adf_mat <- as.matrix(scale(x = adf_mat, center = TRUE, scale = TRUE))
    
  } else if (stand == "max") {
    
    adf_mat <- 
      apply(X = adf_mat, MARGIN = 2, FUN = function(z) { y <- z/max(z)} )
    
  } else { 
    
    stop("this function requires a standardisation method") 
    
  }
  
  
  if (metric == "scaling_mf") {
    
    mf_out <- as.numeric(scale(rowSums(adf_mat), center = TRUE, scale = TRUE))
    
  } else if (metric == "unique_mf") {
    
    mf_out <- apply(adf_mat, MARGIN = 1, function(x) mean(x) ) - apply(adf_mat, MARGIN = 1, function(x) sd(x) )
    
  } else if (metric == "sam_mf") {
    
    mf_out <- apply(adf_mat, MARGIN = 1, function(x) mean(x) )/apply(adf_mat, MARGIN = 1, function(x) sd(x) )
    
  } else { 
    
    stop("this function requires a metric to be chosen") 
    
  }
  
  # add this metric onto the original dataframe
  adf$mf <- mf_out
  
  # change the name of this variable to the name of the metric
  names(adf)[names(adf) == "mf"] <- metric
  
  return(adf)
  
}


# test this function on simulated data
MF_scale_funcs(adf = AvFunc, vars = func.names, 
              stand = "max", 
              metric = "sam_mf")



