
# function to calculate multifunctionality as proposed by Manning et al. (2018, Nature Ecology and Evolution)

# adf is a dataframe with plots in rows, and functions in columns
# vars has to bee a named vector of functions to include which has to correspond to column names
# ind is the index used to determine the optimal number of clusters using NbClust function (if multiple are included, the median of all of them is used)
# met is the method used for clustering. This functions supports all methods in hclust() besides "kmeans"
# thresh is the threshold for function assigngin a zero or one

library(NbClust)

manning_multifunc <- function(adf, vars, 
                              ind = c("mcclain", "cindex", "silhouette", "dunn"), 
                              met = "ward.D2",
                              dis = "euclidean",
                              thresh = 0.5) {
  
  if(! "NbClust" %in% installed.packages()[,1]) stop(
    "this function requires NbClust to be installed and loaded"
  )
  
  if(length(vars) <= 1) stop(
    "this function requires at least two ecosystem functions"
  )
  
  # extract functions from input matrix: adf
  adf_mat <- adf[, vars] 
  
  # transpose the adf matrix
  adf_mat_t <- as.data.frame(t((adf_mat[ , ]))) 
  
  # generate a distance matrix based on the transposed adf matrix
  adf_dist <- dist(adf_mat_t, method = dis)
  
  # from the distance matrix, use NbClust to find optimal number of clusters
  opt_clus <- vector(length = length(ind))
  for(i in 1:length(ind)) {
    
    x <- NbClust::NbClust(diss = adf_dist, 
                          distance = NULL, 
                          min.nc = 2, max.nc = (ncol(adf_mat) - 1), 
                          method = met, index = ind[i])
    
    opt_clus[i] <- x$Best.nc[1]
    
  }
  
  # take the median number of clusters based on the four indices computed from NbClust
  clus_n <- round(x = median(opt_clus), digits = 0)
  
  # use hclust() to cluster the functions
  adf_den <- hclust(func_dist, method = met )
  
  # based on median number of clusters, divide the tree into groups
  adf_groups <- cutree(adf_den, k = clus_n)
  
  # function loadings
  func_loads <- 
    unlist(lapply(adf_groups, FUN = function(x) { ( (1)/length(adf_groups[adf_groups == x]) ) } ))
  
  # standardise functions by maximum and assign based on thresholds
  adf_thres <- 
    apply(X = adf_mat, MARGIN = 2,
          FUN = function(z) { 
            
            y <- z/max(z)
            
            ifelse(y > thresh, 1, 0)
            
          } )
  
  # multiply this by the function loads based on the cluster analysis
  adf_thres <- sweep(adf_thres, MARGIN = 2, func_loads, `*`)
  
  # add the multifunctionality score to the input data
  adf$manning_mf <-  rowSums(adf_thres)/ncol(adf_mat)
  
  return(adf)
  
}






