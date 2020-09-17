
# function to calculate three multifunctionality metric proposed in the literature:
# (1) Pasari et al. (2013) multifunctionality metric
# (2) Dooley (2018) Scaled Average Multifunctionality (SAM)
# (3) Jing et al. (2020) Scaling Mulifunctionality Metric

# adf, is dataframe with plots in rows, and functions in columns
# vars has to bee a named vector of functions to include which has to correspond to column names

adf <- AvFunc

vars <- func.names

stand <- "z_score"


### start the function

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


# Jing et al. (2020, Journal of Plant Ecology): Scaling Multifunctionality Metric

scale(rowSums(adf_mat), center = TRUE, scale = TRUE)[1:5] 

# Pasari et al. (2013, Proceedings of the National Academy of Sciences: Unique Multifunctionality Metric

apply(adf_mat, MARGIN = 1, function(x) mean(x) ) - apply(adf_mat, MARGIN = 1, function(x) sd(x) )

# Dooley (2018): Scaled Average Multifunctionality



