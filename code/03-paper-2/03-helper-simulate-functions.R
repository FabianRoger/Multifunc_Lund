#'
#' @title sim_funcs
#' 
#' @description A function that we can use to simulate a set of linear 
#' functions of the relationships between biodiversity and different
#' ecosystem functions
#' 
#' @param n_func - number of functions to simulate
#' @param n - number of plots to simulate for each function
#' @param lambda - mean of the Poisson distribution to simulate diversity
#' species richness values
#' @param mu_est - mean of the Normal distribution of the slopes of each function
#' @param sd_est - sd of the Normal distribution of the slopes of each function
#' @param error_sd - sd of the Normal distribution for the random error of each function

sim_funcs <- function(n_func = 10, n = 100, 
                      lambda = 10, mu_est = 0.1, sd_est = 0.1,
                      error_sd = 0.5) {
  
  # get the diversity of these different plots
  div <- rpois(n = n, lambda = lambda)
  
  # simulate a set of coefficients for the effect of biodiversity on functioning
  fslope <- rnorm(n = n_func, mean = mu_est, sd_est)
  
  # simulate a set of functions
  fsim_list <-
    
    lapply(fslope, function(x) {
      
      y <- (div*x) + rnorm(n = length(div), mean = 0, sd = error_sd)
      z <- y + abs( min(y) )
      return(z/max(z))
      
    })
  
  # pull the functions into a data.frame
  fsim <- do.call("cbind", fsim_list)
  
  fsim <- list(dplyr::tibble(plot = 1:n, div = div),
               fsim)
  
  return(fsim)
  
}

### END
