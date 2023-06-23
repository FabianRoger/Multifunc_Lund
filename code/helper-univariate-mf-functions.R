
#' @title standardise_functions
#' @description Function to standardise ecosystem functions using several typical standardisation methods
#' 
#' @details Takes a numeric vector and standardises the values using one of five methods:
#' 1. "z_score": (vector - mean(vector))/sd(vector)
#' 2. "z_score_abs": (vector - mean(vector))/sd(vector) + min((vector - mean(vector))/sd(vector))
#' 3. "max_0_1": ( vector - min(vector) )/( max(vector) - min(vector) )
#' 4. "max": vector/max(vector)
#' 5. "max_5_%": vector/max(quantile(vector, 0.95))
#' This function is implemented in all univariate multifunctionality calculations for standardisation
#' 
#' @author James G. Hagan (james_hagan(at)outlook.com)
#' 
#' @param x numeric vector of ecosystem function values
#' @param method method of standardisation ("z_score", "z_score_abs", "max_0_1", "max", "max_5_%", see details above)
#' 
#' @return standardised numeric vector
#' 

standardise_functions <- function(x, method) {
  
  if (method == "z_score") {
    
    y <- scale(x)[,1]
    
  } else if (method == "z_score_abs") {
    
    y <- scale(x)[,1]
    y <- y + abs(min(y))
    
  } else if (method == "max_0_1") {
    
    y <- ( x - min(x) )/( max(x) - min(x) )
    
  } else if (method == "max") {
    
    y <- x/max(x, na.rm = TRUE)
    
  } else if (method == "max_5_%") {
    
    z <- x[ (x != Inf & x != -Inf & !is.na(x)) ]
    max_x <- z[z >= quantile(z, 0.95, na.rm = TRUE)]
    y <- x/mean(max_x)
    
  } else if (method == "none") {
    
    y <- x
    
  } else {
    
    stop("error, choose appropriate standardisation method", call. = FALSE)
    
  }
  
  return(y)
  
}

#' @title cluster_multifunc
#' @description Calculate ecosystem function multifunctionality (Manning et al. 2018) for rows in a data.frame
#' 
#' @details Takes a data frame and ecosystem function names and calculates Manning et al.'s (2018) ecosystem function
#' multifunctionality metric. In this implementation (unlike the original implementation), functions are clustered using
#' an automated clustering procedure. Specifically, we cluster the functions using four algorithms that are implemented
#' in the hclust() function: "mcclain", "cindex", "silhouette", "dunn". The optimal number of clusters is chosen using the
#' "ward.D2" method from the NbClust package. After automated clustering, the method proceeds as described by Manning et al. (2018)
#' 
#' @author James G. Hagan (james_hagan(at)outlook.com)
#' 
#' @param adf dataframe with plots in rows, and functions in columns
#' @param vars character vector with the names of the chosen ecosystem functions corresponding to the column names in adf
#' @param ind clustering methods from the NbClust package, default =  c("mcclain", "cindex", "silhouette", "dunn")
#' @param met method for determining the optimal number of clusters from the NbClust package, only supports "ward.D2" at present
#' @param dis distance metric, default = "euclidean"
#' @param thresh threshold beyond which a function is counted as 1, default = 0.5
#' @param stand_method method of standardisation from the standardise_functions() function (above), default = "max_5_%"
#' 
#' @references 
#' 
#' Manning et al. 2018. Redefining ecosystem multifunctionality. Nature Ecology and Evolution, 2(3): 427-436.
#' 
#' @return returns a vector of ecosystem function multifunctionality (Manning et al. 2018) values for each row
#'
#'

cluster_multifunc <- function(adf, vars, 
                              ind = c("mcclain", "cindex", "silhouette", "dunn"), 
                              met = "ward.D2",
                              dis = "euclidean",
                              thresh = 0.5,
                              stand_method = "max_5_%"
                              ) {
  
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
  
  # scale the transposed matrix
  adf_mat_t <- apply(adf_mat_t, 2, standardise_functions, method = "z_score")
  adf_mat_t <- as.data.frame(adf_mat_t)

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
            
            y <- standardise_functions(x = z, method = stand_method)
            
            ifelse(y > thresh, 1, 0)
            
          } )
  
  # multiply this by the function loads based on the cluster analysis
  adf_thres <- sweep(adf_thres, MARGIN = 2, func_loads, `*`)
  
  # add the multifunctionality score to the input data
  manning_mf <-  rowSums(adf_thres)/clus_n
  
  return(manning_mf)
  
}

#' @title pca_multifunc
#' @description Calculate PCA-based multifunctionality (Meyer et al. 2018) for rows in a data.frame
#' 
#' @details Takes a data frame and ecosystem function names and calculates Meyer et al.'s (2018) PCA-based multifunctionality metric. 
#' The method of standardisation cannot be changed with this method because it needs z-score standardised function data. Otherwise,
#' the method is implemented exactly as per Meyer et al. (2018)
#' 
#' @author James G. Hagan (james_hagan(at)outlook.com) and Fabian Roger
#' 
#' @param adf dataframe with plots in rows, and functions in columns
#' @param vars character vector with the names of the chosen ecosystem functions corresponding to the column names in adf
#' 
#' @references 
#' 
#' Meyer et al. 2018. Biodiversity–multifunctionality relationships depend on identity and number of measured functions. 
#' Nature Ecology and Evolution, 2(1): 44-49.
#' 
#' @return returns a vector of PCA-based multifunctionality (Meyer et al. 2018) values for each row
#'
#'

pca_multifunc <- function(adf, vars){
  
  if(! "vegan" %in% installed.packages()[,1]) stop(
    "this function requires vegan to be installed"
  )
  
  warning("stand_method cannot be altered with this function as the PCA approach requires z_score standardised function data")
  
  adf_mat <- adf[,vars]
  adf_mat <- apply(adf_mat, 2, standardise_functions, method = "z_score")
  adf_mat <- as.data.frame(adf_mat)
  
  pca2<-vegan::rda(adf_mat, scale = FALSE)
  
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

#' @title MF_pasari 
#' @description Calculate Pasari et al.'s (2013) unique multifunctionality metric for rows in a data.frame
#' 
#' @details Takes a data frame, variable names and a standardisation method and calculate the "unique multifunctionality
#' metric" proposed by Pasari et al. (2013). The metric is the average across standardised functions minus the standard deviation
#' among functions. Therefore, the metric decreases when there is a lot of variation among functions.
#' 
#' @author James G. Hagan (james_hagan(at)outlook.com)
#' 
#' @param adf dataframe with plots in rows, and functions in columns
#' @param vars character vector with the names of the chosen ecosystem functions corresponding to the column names in adf
#' @param stand_method method of standardisation from the standardise_functions() function (above), default = "max"
#' 
#' @references 
#' 
#' Pasari, J.R., Levi, T., Zavaleta, E.S. and Tilman, D., 2013. Several scales of biodiversity affect ecosystem 
#' multifunctionality. Proceedings of the National Academy of Sciences, 110(25): 10219-10222.
#' 
#' @return returns a vector of Pasari et al.'s (2013) multifunctionality values for each row
#' 

MF_pasari <- function(adf, vars, stand_method = "max") {
  
  # extract functions from input matrix: adf
  adf_mat <- adf[, vars]
  adf_mat <- apply(adf_mat, 2, standardise_functions, method = stand_method)
  adf_mat <- as.data.frame(adf_mat)
  
  mf_pasari <- apply(adf_mat, MARGIN = 1, function(x) mean(x) ) - apply(adf_mat, MARGIN = 1, function(x) sd(x) ) 
  
  return(mf_pasari)
  
}

#' @title MF_dooley
#' @description Calculate Dooley et al.'s (2018) scaled average multifunctionality metric for rows in a data.frame
#' 
#' @details Takes a data frame, variable names and a standardisation method and calculates the "scaled average multifunctionality
#' metric (SAM metric)" proposed by Dooley et al. (2018). The metric is the mean across standardised functions minus after which
#' the mean across functions is divided by the standard deviation among functions. Therefore, the metric decreases 
#' when there is a lot of variation among functions.
#' 
#' @author James G. Hagan (james_hagan(at)outlook.com)
#' 
#' @param adf dataframe with plots in rows, and functions in columns
#' @param vars character vector with the names of the chosen ecosystem functions corresponding to the column names in adf
#' @param stand_method method of standardisation from the standardise_functions() function (above), default = "max_5_%"
#' 
#' @references 
#' 
#' Dooley, Á., 2018. Modelling Techniques for Biodiversity and Ecosystem Multifunctionality: 
#' Theoretical Development and Application (Doctoral dissertation, National University of 
#' Ireland, Maynooth (Ireland)).
#' 
#' @return returns a vector of Dooley et al.'s (2018) scaled average multifunctionality values for each row
#' 

MF_dooley <- function(adf, vars, stand_method = "max_5_%") {
  
  # extract functions from input matrix: adf
  adf_mat <- adf[, vars] 
  adf_mat <- apply(adf_mat, 2, standardise_functions, method = stand_method)
  adf_mat <- as.data.frame(adf_mat)
  
  mf_dooley <- apply(adf_mat, MARGIN = 1, function(x) mean(x) )/apply(adf_mat, MARGIN = 1, function(x) sd(x) )
  
  return(mf_dooley)
  
}

#' @title MF_jing
#' @description Calculate Jing et al.'s (2020) scaling multifunctionalityu metric for rows in a data.frame
#' 
#' @details Takes a data frame, variable names and a standardisation method and calculates the "scaling multifunctionality
#' metric" proposed by Jing et al. (2020). The metric is calculated by taking the mean among z-score-standardised functions and then
#' standardising the mean multifunctionality values using the z-score standardisation method. This method can only be implemented
#' using the z-score standardisation because the method only makes sense using this standardisation.
#' 
#' @author James G. Hagan (james_hagan(at)outlook.com)
#' 
#' @param adf dataframe with plots in rows, and functions in columns
#' @param vars character vector with the names of the chosen ecosystem functions corresponding to the column names in adf
#' 
#' @references 
#' 
#' Jing, X., Prager, C.M., Classen, A.T., Maestre, F.T., He, J.S. and Sanders, N.J., 2020. Variation in the methods 
#' leads to variation in the interpretation of biodiversity–ecosystem multifunctionality relationships. 
#' Journal of Plant Ecology, 13(4): 431-441.
#' 
#' @return returns a vector of Jing et al.'s (2018) scaling multifunctionality metric values for each row
#' 

MF_jing <- function(adf, vars) {
  
  # extract functions from input matrix: adf
  adf_mat <- adf[, vars] 
  adf_mat <- apply(adf_mat, 2, standardise_functions, method = "z_score")
  adf_mat <- as.data.frame(adf_mat)
  
  warning("stand_method cannot be altered with this function as Jing et al.'s approach requires z_score standardised function data")
  
  mf_jing <- standardise_functions(x = rowSums(adf_mat), method = "z_score")
  
  return(mf_jing) 
  
}


# summing approach

# function to calculate summed multifunctionality

# adf, is dataframe with plots in rows, and functions in columns
# vars has to bee a named vector of functions to include which has to correspond to column names
MF_sum <- function(adf, vars, stand_method = "z_score_abs") {
  
  adf_mat <- adf[, vars]
  adf_mat <- apply(adf_mat, 2, standardise_functions, method = stand_method)
  adf_mat <- as.data.frame(adf_mat)
  
  mf_sum <- rowSums(adf_mat)
  
  return(mf_sum)
  
}

# average approach

# function to calculate average multifunctionality

# adf, is dataframe with plots in rows, and functions in columns
# vars has to bee a named vector of functions to include which has to correspond to column names
# stand_method = method used to standardise the data ("none", "z_score", "z_score_abs", "max", "max_0_1", "max_5_%")

MF_av <- function(adf, vars, stand_method = "z_score") {
  
  adf_mat <- adf[, vars]
  adf_mat <- apply(adf_mat, 2, standardise_functions, method = stand_method)
  adf_mat <- as.data.frame(adf_mat)
  
  mf_av <- rowSums(adf_mat)/length(vars)
  
  return(mf_av)
  
}

# function to calculate standard deviation among functions

# adf, is dataframe with plots in rows, and functions in columns
# vars has to bee a named vector of functions to include which has to correspond to column names
# stand_method = method used to standardise the data ("none", "z_score", "z_score_abs", "max", "max_0_1", "max_5_%")

MF_sd <- function(adf, vars, stand_method = "z_score_abs") {
  
  adf_mat <- adf[, vars]
  adf_mat <- apply(adf_mat, 2, standardise_functions, method = stand_method)
  adf_mat <- as.data.frame(adf_mat)
  
  mf_sd <- apply(adf_mat, 1, sd, na.rm = TRUE)
  
  return(mf_sd)
  
}

# function to calculate coefficient of variation among functions

# adf, is dataframe with plots in rows, and functions in columns
# vars has to bee a named vector of functions to include which has to correspond to column names
# stand_method = method used to standardise the data ("none", "z_score", "z_score_abs", "max", "max_0_1", "max_5_%")

MF_cv <- function(adf, vars, stand_method = "z_score_abs") {
  
  adf_mat <- adf[, vars]
  adf_mat <- apply(adf_mat, 2, standardise_functions, method = stand_method)
  adf_mat <- as.data.frame(adf_mat)
  
  mf_cv <- apply(adf_mat, 1, function(x) { (sd(x)/mean(x)) })
  
  return(mf_cv)
  
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
  
  return(funcMaxed)
  
}


# Simpson's diversity index: Raudsepp-Hearne et al. (2009)

# adf, is dataframe with plots in rows, and functions in columns
# vars has to bee a named vector of functions to include which has to correspond to column names
# stand_method = method used to standardise the data ("none", "z_score", "z_score_abs", "max", "max_0_1", "max_5_%")

MF_simpsons_div <- function(adf, vars, stand_method = "max") {
  
  if(! "vegan" %in% installed.packages()[,1]) stop(
    "this function requires vegan to be installed"
  )
  
  adf_mat <- adf[, vars]
  adf_mat <- apply(adf_mat, 2, standardise_functions, method = stand_method)
  adf_mat <- as.data.frame(adf_mat)
  
  mf_simp <- vegan::diversity(x = adf_mat, index = "simpson", MARGIN = 1)
  
  z <- 
    sapply(adf_mat, function(x) {
    
    ifelse(sum(x) == 0, 1, 0)
    
  })
  
  if (sum(z) == length(vars)) {
    
    warning(x = "all functions are zero which leads to maximum Simpson diversity")
    
    return(mf_simp)
    
  } else {
    
    return(mf_simp)
    
  }
  
}

# Shannon diversity index

# adf, is dataframe with plots in rows, and functions in columns
# vars has to bee a named vector of functions to include which has to correspond to column names
# stand_method = method used to standardise the data ("none", "z_score", "z_score_abs", "max", "max_0_1", "max_5_%")

MF_shannon_div <- function(adf, vars, stand_method = "max") {
  
  if(! "vegan" %in% installed.packages()[,1]) stop(
    "this function requires vegan to be installed"
  )
  
  adf_mat <- adf[, vars]
  adf_mat <- apply(adf_mat, 2, standardise_functions, method = stand_method)
  adf_mat <- as.data.frame(adf_mat)
  
  mf_shan <- vegan::diversity(x = adf_mat, index = "shannon", MARGIN = 1)
  
  return(mf_shan)
  
}


# MESLI approach Rodríguez-Loinaz et al. (2015)

# adf, is dataframe with plots in rows, and functions in columns
# vars has to bee a named vector of functions to include which has to correspond to column names

MF_mesli <- function(adf, vars, stand_method = "max_0_1") {
  
  adf_mat <- adf[, vars]
  adf_mat <- apply(adf_mat, 2, standardise_functions, method = stand_method)
  adf_mat <- as.data.frame(adf_mat)
  
  mf_mesli <- rowSums(adf_mat)
  
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
    
    a <- apply(X = adf_mat, MARGIN = 2, function(x) { mean(x[x <= quantile(x, probs = A_quant)], na.rm = TRUE)   })
    
  } else {
    
    stop("this function requires correct specification of A values, see Documentation")
    
  }
  
  # set up the B values
  if ( ((is.vector(B) == TRUE) & (length(B) == length(vars)))  ) {
    
    b <- B
    
  } else if (B == "max") {
    
    b <- apply(X = adf_mat, MARGIN = 2, function(x) { mean(x[x >= quantile(x, probs = B_quant)], na.rm = TRUE) })
    
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
  
  # apply the geometric mean
  d_mf <- 
    apply(d_out, MARGIN = 1, function(x) {
      
      exp( mean(log(x)) )
      
    })
  
  d_mf <- (d_mf)^(1/sum(ws))
  
  return(d_mf)
  
}

### END
