
# functions to calculate multifunctionality

# load the relevant libraries
library(vegan)
library(NbClust)


# key arguments for all functions:
# adf, is dataframe with plots in rows, and functions in columns
# vars has to be a named vector


# Hill approach (Roger et al. unpublished)

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


# Manning et al. (2018) approach

# function to calculate multifunctionality as proposed by Manning et al. (2018, Nature Ecology and Evolution)

# adf is a dataframe with plots in rows, and functions in columns
# vars has to bee a named vector of functions to include which has to correspond to column names
# ind is the index used to determine the optimal number of clusters using NbClust function (if multiple are included, the median of all of them is used)
# met is the method used for clustering. This functions supports all methods in hclust() besides "kmeans"
# thresh is the threshold for function assigngin a zero or one

# this doesn't work with few functions

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
                          min.nc = 2, max.nc = (ncol(adf_mat)-1), 
                          method = met, index = ind[i])
    
    opt_clus[i] <- x$Best.nc[1]
    
  }
  
  # take the median number of clusters based on the four indices computed from NbClust
  clus_n <- round(x = median(opt_clus), digits = 0)
  
  # use hclust() to cluster the functions
  adf_den <- hclust(adf_dist, method = met )
  
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


# Meyer et al. (2017) approach

# function to calculate the pca based multifunctionality index as suggested by 
# Meyer et al. (2017) Biodiversity-multifunctionality relationships depend on identity and number of measured functions #
# Nature Ecology & Evolution   

# adf, is dataframe with plots in rows, and functions in columns
# vars has to bee a named vector
# standardise is a logical value (TRUE or FALSE)

pca_multifunc <- function(adf, vars, standardise = FALSE){
  
  if(! "vegan" %in% installed.packages()[,1]) stop(
    "this function requires vegan to be installed"
  )
  
  adf_mat <- adf[,vars]
  pca2<-vegan::rda(adf_mat, scale = standardise)
  
  temp2 <- vegan::scores(pca2, choices=1:length(vars), display=c("sites")) 
  eig<-summary(pca2)$cont$importance[1,]
  for(i in 1:length(eig)) temp2[,i] <- temp2[,i] * eig[i]
  Index.wt <- rowSums(temp2)
  adf$multifunc_pca_ind <-  Index.wt
  return(adf)
  
}


# Pasari et al. (2013) approach

# function to calculate three multifunctionality metric proposed in the literature:
# Pasari et al. (2013) multifunctionality metric

# adf, is dataframe with plots in rows, and functions in columns
# vars has to bee a named vector of functions to include which has to correspond to column names
# stand is the method of standardisation. three methods are supported: (1) no standardisation ("none") (2) z-score standardisation ("z_score") and (3) by the maximum ("max")

MF_pasari <- function(adf, vars, stand = "none") {
  
  # extract functions from input matrix: adf
  adf_mat <- adf[, vars] 
  
  # standardise the functions using: (1) the z-score ("z_score") or (2) by the maximum ("max")
  
  if (stand == "none") {
    
  } else if (stand == "z_score") {
    
    adf_mat <- as.matrix(scale(x = adf_mat, center = TRUE, scale = TRUE))
    
  } else if (stand == "max") {
    
    adf_mat <- 
      apply(X = adf_mat, MARGIN = 2, FUN = function(z) { z/max(z)} )
    
  } else { 
    
    stop("this function requires a standardisation method") 
    
  }
  
  adf$mf_pasari <- apply(adf_mat, MARGIN = 1, function(x) mean(x) ) - apply(adf_mat, MARGIN = 1, function(x) sd(x) ) 
  
  return(adf)
  
}


# Dooley et al. (2018) approach

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
      apply(X = adf_mat, MARGIN = 2, FUN = function(z) { z/max(z)} )
    
  } else { 
    
    stop("this function requires a standardisation method") 
    
  }
  
  adf$mf_dooley <- apply(adf_mat, MARGIN = 1, function(x) mean(x) )/apply(adf_mat, MARGIN = 1, function(x) sd(x) )
  
  return(adf)
  
}


# Jing et al. (2020) approach

# function to calculate three multifunctionality metric proposed in the literature:
# Jing et al. (2020) Scaling Mulifunctionality Metric

# adf, is dataframe with plots in rows, and functions in columns
# vars has to bee a named vector of functions to include which has to correspond to column names
# stand is the method of standardisation. three methods are supported: (1) no standardisation ("none") (2) z-score standardisation ("z_score") and (3) by the maximum ("max")

MF_jing <- function(adf, vars, stand = "none") {
  
  # extract functions from input matrix: adf
  adf_mat <- adf[, vars] 
  
  # standardise the functions using: (1) the z-score ("z_score") or (2) by the maximum ("max")
  
  if (stand == "none") {
    
  } else if (stand == "z_score") {
    
    adf_mat <- as.matrix(scale(x = adf_mat, center = TRUE, scale = TRUE))
    
  } else if (stand == "max") {
    
    adf_mat <- 
      apply(X = adf_mat, MARGIN = 2, FUN = function(z) { z/max(z)} )
    
  } else { 
    
    stop("this function requires a standardisation method") 
    
  }
  
  adf$mf_jing <- as.numeric(scale(rowSums(adf_mat), center = TRUE, scale = TRUE))
  
  return(adf) 
  
}


# test these functions

# generate the simulated data
set.seed(777)
specnum <- 10
funcnum <- 10
distribution = "runif"
FuncMat <- FunctionValue(specnum,funcnum, distribution, min = 0.1, max = 0.9)
func.names <- as.character( unique( FuncMat$Functions))
spec.names <- as.character( unique( FuncMat$Species))
maxrep <- 100 
SpecMat <- SpeciesMatrix(specnum = specnum, maxrep = maxrep)
method = "av"
compfunc <- func.names[1:3]
AvFunc <- AverageFunction(SpecMat, FuncMat,
                          method = method, 
                          compfunc = compfunc)
set.seed(563)
errM <- matrix(rnorm(n = nrow(AvFunc)*funcnum, mean = 0, sd = 0.01), ncol = funcnum)
AvFunc[,func.names] <- AvFunc[,func.names] + errM

# simulated dataset
AvFunc


# Hill approach
hill_multifunc(adf = AvFunc, vars = func.names, scale = 1, HILL = TRUE)

# Manning approach
manning_multifunc(adf = AvFunc, vars = func.names, 
                  ind = c("mcclain", "cindex", "silhouette", "dunn"), 
                  met = "ward.D2",
                  dis = "euclidean",
                  thresh = 0.5)

# Meyer approach
pca_multifunc(adf = AvFunc, vars = func.names, standardise = FALSE)

# Pasari approach 
MF_pasari(adf = AvFunc, vars = func.names, stand = "none")

# Dooley approach
MF_dooley(adf = AvFunc, vars = func.names, stand = "max") 

# Jing approach
MF_jing(adf = AvFunc, vars = func.names, stand = "max") 










