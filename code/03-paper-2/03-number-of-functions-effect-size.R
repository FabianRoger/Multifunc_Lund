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
source("code/03-paper-2/helper-simulate-functions.R")
source("code/03-paper-2/helper-slope-estimate-slope-n-functions.R")
source("code/helper-plotting-theme.R")
source("code/helper-univariate-mf-functions.R")

# set the number of simulations to do
metrics <- c("ave", "inv-Simpson", "thresh_30", "thresh_70", "ENFQ1", "ENFQ2")
title <- c("Average MF", "Inv-Simpson MF", "Thresh 30% MF",  "Thresh 70% MF",
           "ENF-Q1 MF", "ENF-Q2 MF")
xlabs <- c(NA, NA, NA, NA, "Number of functions", "Number of functions")
ylabs <- c("Slope est. (+-SE)", NA, "Slope est. (+-SE)", NA, "Slope est. (+-SE)")

# set the number of simulations
n_sim <- 100

# random functions or linked functions
random <- FALSE

if(random) {
  
  # simulate the datasets
  fsim_list <- 
    lapply(1:n_sim, function(x) {
      
      fsim <- sim_funcs(n_func = 9, n = 100, 
                        lambda = 10, 
                        mu_est = 0.1, sd_est = 0.1, 
                        error_sd = 0.5)
      
      # convert to adf_dat
      fsim[[2]] <- dplyr::as_tibble(fsim[[2]])
      names(fsim[[2]]) <- paste0("F", 1:ncol(fsim[[2]]))
      
      return(fsim)
      
    })
  
} else {
  
  # random functions
  fsim_list <- 
    lapply(1:n_sim, function(x) {
      
      fsim <- sim_funcs(n_func = 3, n = 100, 
                        lambda = 10, 
                        mu_est = 0.1, sd_est = 0.1, 
                        error_sd = 0.5)
      
      fsim2 <- apply(fsim[[2]], 2, function(x) { 
        
        y <- (x*rnorm(n = 1, mean = 0.1, 0.1)) + rnorm(n = length(x), 0, 0.01)
        z <- y + abs( min(y) )
        z/max(z)
        
      } )
      
      # convert to adf_dat
      fsim[[2]] <- dplyr::as_tibble( cbind( fsim[[2]], fsim2) )
      names(fsim[[2]]) <- paste0("F", 1:ncol(fsim[[2]]))
      
      return(fsim)
      
    })
  
}

plot_list <- vector("list", length = length(metrics))
for(i in 1:length(metrics)){
  
  # plot an example figure
  
  # get a random dataset
  df_ex <- fsim_list[[sample(x = 1:n_sim, size = 1)]]
  
  # convert to adf_dat
  ex_nfunc <- 
    n_func_est(adf = df_ex[[2]], vars = names(df_ex[[2]]), 
               div = df_ex[[1]]$div, 
               metric = metrics[i])
  
  # add a est id variable
  ex_nfunc$est_id <- with(ex_nfunc, paste(n_func, id, sep = "_"))
  
  ex_nfunc_sum <- 
    ex_nfunc |>
    dplyr::group_by(n_func) |>
    dplyr::summarise(slope_est_m = mean(slope_est), .groups = "drop")
  
  p1 <- 
    ggplot(data = ex_nfunc) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "black") +
    geom_point(mapping = aes(x = n_func, y = slope_est, group = est_id), 
               shape = 1, alpha = 0.2, position = position_dodge(width = 0.2)) +
    geom_errorbar(mapping = aes(x = n_func, 
                                ymin = slope_est - slope_SE, ymax = slope_est + slope_SE,
                                group = est_id),
                  width = 0, alpha = 0.2, position = position_dodge(width = 0.2)) +
    geom_line(data = ex_nfunc_sum, 
              mapping = aes(x = n_func, y = slope_est_m), colour = "red") +
    ylab(if(is.na(ylabs[i])){NULL}else{ylabs[i]}) +
    xlab(if(is.na(xlabs[i])){NULL}else{xlabs[i]}) +
    scale_x_continuous(breaks = c(2:9)) +
    ggtitle(title[i]) +
    theme_meta() +
    theme(plot.title = element_text(size = 11),
          axis.text.x = element_text(size=10))
  
  # perform the analysis for all simulated datasets
  nfunc_plot <- 
    
    lapply(fsim_list, function(x) {
      
      # get the function data
      adf_dat <- x[[2]]  
      
      # estimate the slope between div and EMF with different numbers of functions
      df_nfunc <- 
        n_func_est(adf = adf_dat, vars = names(adf_dat), 
                   div = x[[1]]$div, 
                   metric = metrics[i])
      
      df_nfunc <- 
        df_nfunc |>
        dplyr::group_by(n_func) |>
        dplyr::summarise(slope_est_m = mean(slope_est))
      
      return(df_nfunc)
      
    } )
  
  # bind into a data.frame
  nfunc_plot <- dplyr::bind_rows(nfunc_plot, .id = "id")
  
  # plot the different slopes
  p2 <- 
    ggplot(data = nfunc_plot) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "black") +
    geom_line(mapping = aes(x = n_func, y = slope_est_m, group = id), 
              colour = "red", linewidth = 0.25, alpha = 0.5) +
    ylab(NULL) +
    xlab(if(is.na(xlabs[i])){NULL}else{xlabs[i]}) +
    scale_x_continuous(breaks = c(2:9)) +
    ggtitle(" ") +
    theme_meta() +
    theme(plot.title = element_text(size = 11),
          axis.text.x = element_text(size=10))
  
  # join the two plots
  p12 <- cowplot::plot_grid(p1, p2, rel_widths = c(1, 1.15),
                            align = "h")
  
  # add the plot to a list
  plot_list[[i]] <- p12
  
}

# combine into a single plot
p <- 
  cowplot::plot_grid(plotlist = plot_list,
                     nrow = 3, ncol = 2, labels = letters[1:6],
                     label_size = 11, label_fontface = "plain",
                     rel_heights = c(1, 1, 1.075),
                     rel_widths = c(1.085, 1, 1),
                     align = "h")
plot(p)

if(random) {
  
  ggsave("figures-paper-2/fig_2.svg", p,
         units = "cm", height = 18, width = 23)
  
} else {
  
  ggsave("figures-paper-2/fig_S1.svg", p,
         units = "cm", height = 18, width = 23)
  
}

### END
