
# function to calculate the effective number of functions as proposed by Roger, Bagchi and Byrnes (unpublished!) 

# adf, is dataframe with plots in rows, and functions in columns
# vars has to bee a named vector of functions to include which has to correspond to column names
# scales is the order of diversity that should be calculated. default = 1

hill_multifunc <- function(adf, vars, scale = 1, HILL = TRUE){
  
  if(! "vegan" %in% installed.packages()[,1]) stop(
    "this function requires vegan to be installed"
  )
  
  if(length(scale) > 1) stop(
    "this function requires to choose a single order of diversity"
  )
  
  adf_mat <- adf[,vars]
  effN <- vegan::renyi(adf_mat, scales = scale, hill = TRUE)
  if(!HILL) effN <- effN/length(vars)
  meanFunc <- rowMeans(adf_mat) 
  
  adf$multifunc_effN <-  effN*meanFunc
  return(adf)
  
}
