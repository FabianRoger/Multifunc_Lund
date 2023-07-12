#'
#' @title Test whether slope changes with number of functions
#' 
#' @description Here we simulate relationships between biodiversity and
#' ecosystem multifunctionality (EMF) and then we test whether the relationship
#' between biodiversity and EMF changes by considering different numbers
#' of functions.
#'

# load relevant libraries
library(dplyr)
library(tidyr)
library(ggplot2)

# load relevant functions
source("code/03-paper-2/02-helper-simulate-functions.R")
source("code/helper-plotting-theme.R")
source("code/helper-univariate-mf-functions.R")

# simulate a set of functions
fsim <- sim_funcs(n_func = 3, n = 100, 
                  lambda = 10, 
                  mu_est = 0.25, sd_est = 0.5, 
                  error_sd = 0.5)

# bind the simulation data for efficient plotting
fsim_plot <- dplyr::bind_cols( fsim[[1]], fsim[[2]])

# rename the columns
names(fsim_plot) <- c("plot", "div", paste0("F", 1:ncol(fsim[[2]])) )

# pull into the long format
fsim_plot <- 
  fsim_plot |>
  tidyr::pivot_longer(cols = contains("F"),
                      names_to = "Function", 
                      values_to = "Value")

p1 <- 
  ggplot(data = fsim_plot,
         mapping = aes(x = div, y = Value, group = Function)) +
  geom_smooth(method = "lm", size = 0.5, colour = "grey", se = FALSE) +
  ylab("Function value (0-1)") +
  xlab("Species richness") +
  theme_test() +
  theme(axis.text = element_text(colour = "black"))
plot(p1)

adf_dat <- dplyr::as_tibble(fsim[[2]])
names(adf_dat) <- paste0("F", 1:ncol(fsim[[2]]))

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
          
          multifunc::eff_num_func(dat = adf, vars = func_vec)
          
        } else if(metric == "Shannon") {
            
          MF_shannon(adf = adf, vars = func_vec, stand_method = "none")
          
        } else {
            
          stop("choose appropriate EMF metric")
          
          }
      
      # fit a linear model
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

n_func_est(adf = adf_dat, vars = names(adf_dat), div = fsim[[1]]$div, metric = "thresh_30")


# add a few non-linear functions




