
# functions to calculate multifunctionality

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
  
  multifunc_effN <-  effN*meanFunc
  return( as.numeric(multifunc_effN) )
  
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
  manning_mf <-  rowSums(adf_thres)/clus_n
  
  return(manning_mf)
  
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
  
  pc_load <- prcomp(x = adf_mat)
  pc_load <- pc_load$rotation %*% diag(pc_load$sdev) 
  
  inv <- 
    apply(X = pc_load, MARGIN = 2, function(z) { 
      
      ifelse(z[(max(abs(z)) == z | max(abs(z)) == -z)] > 0, 1, -1)
      
    })
  
  temp2 <- vegan::scores(pca2, choices=1:length(pca2$CA$eig), display=c("sites")) 
  eig<-summary(pca2)$cont$importance[1,]
  for(i in 1:length(eig)) temp2[,i] <- temp2[,i] * eig[i]
  
  temp2 <- sweep(temp2, MARGIN = 2, inv, `*`)
  
  Index.wt <- rowSums(temp2)
  multifunc_pca_ind <-  Index.wt
  return(multifunc_pca_ind)
  
}


# Pasari et al. (2013) approach

# function to calculate Pasari et al.'s (2013) multifunctionality metric

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
  
  mf_pasari <- apply(adf_mat, MARGIN = 1, function(x) mean(x) ) - apply(adf_mat, MARGIN = 1, function(x) sd(x) ) 
  
  return(mf_pasari)
  
}


# Dooley et al. (2018) approach

# function to calculate Dooley's (2018) Scaled Average Multifunctionality (SAM)

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
  
  mf_dooley <- apply(adf_mat, MARGIN = 1, function(x) mean(x) )/apply(adf_mat, MARGIN = 1, function(x) sd(x) )
  
  return(mf_dooley)
  
}


# Jing et al. (2020) approach

# function to calculate Jing et al.'s (2020) Scaling Mulifunctionality Metric

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
  
  mf_jing <- as.numeric(scale(rowSums(adf_mat), center = TRUE, scale = TRUE))
  
  return(mf_jing) 
  
}


# summing approach

# function to calculate summed multifunctionality

# adf, is dataframe with plots in rows, and functions in columns
# vars has to bee a named vector of functions to include which has to correspond to column names
MF_sum <- function(adf, vars) {
  
  adf_mat <- adf[, vars]
  
  mf_sum <- rowSums(adf_mat)
  
  return(mf_sum)
  
}

# average approach

# function to calculate average multifunctionality

# adf, is dataframe with plots in rows, and functions in columns
# vars has to bee a named vector of functions to include which has to correspond to column names
MF_av <- function(adf, vars) {
  
  adf_mat <- adf[, vars]
  
  mf_av <- rowSums(adf_mat)/length(vars)
  
  return(mf_av)
  
}


# single threshold approach

# The function below taken from the multifunc package and written by Jarret Byrnes (https://github.com/jebyrnes/multifunc)

# They are included here in a modified form as 
# 1) the package is not available on CRAN for the moment
# 2) the package loads plyr which conflicts with dplyr if loaded afterwards

# note that I re-wrote the functions below in order to not rely on plyr

single_threshold_mf <- function(adf, vars = NA, thresh = 0.7, prepend = "Diversity", maxN = 1){
  
  if(is.na(vars)[1]) stop("You need to specify some response variable names")
  
  getMaxValue <- function(x){
    l <- length(x)    
    mean( sort(x, na.last=F)[l:(l-maxN+1)], na.rm=T)
  }
  
  funcMaxed <- rowSums(apply(adf[,which(names(adf)%in%vars)], 2, function(x) x >= thresh*getMaxValue(x)))
  
  funcMaxed
  
}


# Simpson's diversity index: Raudsepp-Hearne et al. (2009)

# adf, is dataframe with plots in rows, and functions in columns
# vars has to bee a named vector of functions to include which has to correspond to column names

MF_simpsons_div <- function(adf, vars) {
  
  if(! "vegan" %in% installed.packages()[,1]) stop(
    "this function requires vegan to be installed"
  )
  
  adf_mat <- adf[, vars]
  
  mf_simp <- vegan::diversity(x = adf_mat, index = "simpson", MARGIN = 1)
  
  z <- 
    sapply(adf_mat, function(x) {
    
    ifelse(sum(x) == 0, 1, 0)
    
  })
  
  if (sum(z) == length(vars)) {
    
    warning(x = "all functions are zero which leads to maximum Simpson diversity")
    
    mf_simp
    
  } else {
    
    mf_simp
    
  }
  
}


# MESLI approach Rodríguez-Loinaz et al. (2015)

# adf, is dataframe with plots in rows, and functions in columns
# vars has to bee a named vector of functions to include which has to correspond to column names

MF_mesli <- function(adf, vars) {
  
  adf_mat <- adf[, vars]
  
  y <- 
    apply(adf_mat, MARGIN = 2, FUN = function(x) {
      
      (x - min(x))/(max(x) - min(x))
      
    } )
  
  mf_mesli <- rowSums(y)
  
  return(mf_mesli)

}


# Slade et al. (2017) desirability function approach

# function to calculate multifunctional desirability sensu Slade et al. (2017)

# adf, is dataframe with plots in rows, and functions in columns
# vars has to bee a named vector of functions to include which has to correspond to column names
# A is the lower limit of function deemed acceptable
# B is the upper limit of the function
# s is a shape parameter describing how multifunctional desirability changes between A and B
# when s = 1 (default), multifunctional desirability increases linearly from A to B
# when s < 1,  multifunctional desirability is higher at lower levels of the response
# when s > 1 it is more important to be as close as possible to higher multifunctional desirability
# weights are direct weightings applied to each function

# A, B, s and weights must be defined for each function and must be inputted as vectors (one number for each function)
# if vectors are not inputted for these quantities the defaults are:
# A: for each function, the mean of all values below a certain quantile is used (argument = A_quant)
# B: for each function, the mean of all values above a certain quantile is used (argument = B_quant)
# s: by default, s is 1 which means there is a linear relationship
# weights: equal weights (i.e. weights = 1 for all functions) are assigned by default

MF_slade <- function(adf, vars, 
                     A = "min", 
                     A_quant = 0.10,
                     B = "max", 
                     B_quant = 0.90,
                     s = "linear",
                     weights = "equal") {
  
  
  adf_mat <- adf[, vars]
  
  # set up the A values
  if ( ((is.vector(A) == TRUE) & (length(A) == length(vars)))  ) {
    
    a <- A
    
  } else if (A == "min") {
    
    a <- apply(X = adf_mat, MARGIN = 2, function(x) { mean(x[x < quantile(x, probs = A_quant)], na.rm = TRUE)   })
    
  } else {
    
    stop("this function requires correct specification of A values, see Documentation")
    
  }
  
  # set up the B values
  if ( ((is.vector(B) == TRUE) & (length(B) == length(vars)))  ) {
    
    b <- B
    
  } else if (B == "max") {
    
    b <- apply(X = adf_mat, MARGIN = 2, function(x) { mean(x[x > quantile(x, probs = B_quant)], na.rm = TRUE) })
    
  } else {
    
    stop("this function requires correct specification of B values, see Documentation")
    
  }
  
  
  # set up the s-values
  if ( ((is.vector(s) == TRUE) & (length(s) == length(vars)))  ) {
    
    S <- s
    
  } else if (s == "linear") {
    
    S <- rep(x = 1, times = length(vars))
    
  } else {
    
    stop("this function requires correct specification of s values, see Documentation")
    
  }
  
  d_out <- vector("list", length = length(vars))
  
  for (i in 1:length(vars)) {
    
    p <- (adf_mat[[i]])
    
    y <- ( ( p - a[[i]] )/( b[[i]] - a[[i]] ) )^S[[i]]
    
    z <- mapply(function(d, e) {ifelse(d < a[[i]], 0, e)}, p, y  )
    
    w <- mapply(function(f, g) {ifelse(f > b[[i]], 1, g)}, p, z  )
    
    d_out[[i]] <- w
    
  }
  
  d_out <- do.call("cbind", d_out)
  
  
  # apply weights to each function
  if ( ((is.vector(weights) == TRUE) & (length(weights) == length(vars)))  ) {
    
    ws <- weights
    
  } else if (weights == "equal") {
    
    ws <- rep(x = 1, times = length(vars))
    
  } else {
    
    stop("this function requires correct specification of function weights, see Documentation")
    
  }
  
  # apply the weights by adding it as an exponent
  d_out <- sweep(d_out, MARGIN = 2, ws, `^`)
  
  d_mf <- 
    apply(d_out, MARGIN = 1, function(x) {
      
      exp( mean(log(x)) )
      
    })
  
  d_mf <- (d_mf)^(1/sum(ws))
  
  return(d_mf)
  
}


# define a function to calculate different multifunctionality metrics
# mf.functions is a vector of the names of the function to call
# mf.names are the names to assign to the output of each function call
# add.args is a list of additional arguments to pass to each function called
# if no additional arguments are required then an NA must be specified
multifunc_calculator <- 
  function(adf,
           vars,
           mf.functions = c("MF_sum", "MF_av", "MF_pasari", "single_threshold_mf", "single_threshold_mf"),
           mf.names = c("sum_MF", "ave._MF", "Pasari_MF", "thresh.30_MF", "thresh.70_MF"),
           add.args = list(NA, NA, NA, c(thresh = 0.3), c(thresh = 0.7))) {
    
    # perform error checks
    if( !any(is.na(add.args)) ) {
      print("warning: no NAs in add.args list, is this correct?")
    }
    
    if( length(mf.functions) != length(mf.names) ) {
      stop("error, there must be a name for each multifunctionality metric")
    }
    
    if(length(vars) < 2) {
      stop("error, there must be two or more functions")
    }
    
    # add function names to the additional argument vector
    names(add.args) <- mf.names
    
    # for each function, create an input list with or without additional arguments
    input.list <- vector("list")
    for(j in 1:length(add.args)) {
      
      if (!is.na(add.args[[j]])){
        
        input.list[[j]] <- list(adf, vars)
        
        for (i in 1:length(add.args[[j]])){
          
          input.list[[j]] <- c(input.list[[j]], add.args[[j]][i])
          
        }
        
      } else(
        
        input.list[[j]] <- list(adf, vars)
        
      )
      
    }
    
    # assign each function name the multifunctionality value
    for (k in 1:length(mf.names)) {
      assign(mf.names[k], do.call(mf.functions[k], input.list[[k]]))
    }
    
    # pull the metrics into a data.frame
    mf.df <- sapply(mf.names, function(x){get(x)})
    mf.df <- data.frame(mf.df)
    
    # bind the metrics into a data.frame
    return( data.frame(cbind(adf, mf.df)) )
    
  }


# function to standardise functions and translate them by minimum absolute value
standardise <- function(x) {
  mean.x <- mean(x)
  sd.x <- sd(x)
  x.standardised <- ((x-mean.x)/sd.x)
  x.standardised.positive <- ( x.standardised + abs(min(x.standardised)) )
  return(x.standardised.positive)
}


### END
