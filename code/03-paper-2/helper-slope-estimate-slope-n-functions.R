#'
#' @title n_func_est
#' 
#' @description Function that takes a data.frame or tibble of functions and
#' a vector of biodiversity values and calculates the slope between biodiversity
#' and EMF using different metrics of multifunctionality and considering
#' different numbers of functions
#' 
#' @param adf - dataframe with plots in rows, and functions in columns
#' @param vars - character vector with the names of the chosen ecosystem functions 
#' corresponding to the column names in adf
#' @param div - vector of biodiversity values
#' @param metric - EMF metric: "ave", "thresh_30", "thresh_70", "ENFQ1", "ENFQ2", "inv-Simpson"
#' 

n_func_est <- function(adf, vars, div, metric = "ave") {
  
  func_comb <- vector("list", length = (length(vars)-1) )
  for (i in 2:length(vars)){
    func_comb[[i-1]] <- combn(x = vars, m = i)
  }
  
  list_out <- vector("list", length = length(func_comb))
  for(j in 1:length(func_comb)) {
    
    # get the combinations with j+1 functions
    func_x <- func_comb[[j]]
    
    # get the number of combinations of j+1 functions
    combs <- ncol(func_x)
    
    # vector of output slopes
    comb_B <- vector(length = combs)
    comb_SE <- vector(length = combs)
    for(k in 1:combs) {
      
      # get the set of functions
      func_vec <- func_x[,k]
      
      # calculate EMF
      EMF <- 
        if(metric == "ave") {
          
          MF_av(adf = adf, vars = func_vec, stand_method = "none")
          
        } else if(metric == "thresh_30") {
          
          MF_thresh(adf = adf, vars = func_vec, thresh = 0.3)/length(func_vec)
          
        } else if(metric == "thresh_70") {
          
          MF_thresh(adf = adf, vars = func_vec, thresh = 0.7)/length(func_vec)
          
        } else if(metric == "ENFQ1") {
          
          multifunc::getMF_eff(data = adf, vars = func_vec, q = 1,
                               standardized = TRUE) 
          
        } else if(metric == "ENFQ2") {
          
          multifunc::getMF_eff(data = adf, vars = func_vec, q = 2,
                               standardized = TRUE) 
          
        } else if(metric == "inv-Simpson") {
          
          MF_inv_simpson(adf = adf, vars = func_vec, stand_method = "none")
          
        } else {
          
          stop("choose appropriate EMF metric")
          
        }
      
      # fit a linear model on unstandardised data
      x <- lm(EMF ~ div)
      y <- summary(x)
      comb_B[k] <- y$coefficients[2,][["Estimate"]]
      comb_SE[k] <- y$coefficients[2,][["Std. Error"]]
      
    }
    
    list_out[[j]] <- dplyr::tibble(n_func = j+1,
                                   id = 1:length(comb_B),
                                   slope_est = comb_B,
                                   slope_SE = comb_SE)
    
  }
  
  # bind into a data.frame
  df_est <- dplyr::bind_rows(list_out)
  
  return(df_est)
  
}

### END
